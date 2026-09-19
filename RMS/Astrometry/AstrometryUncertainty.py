""" Uncertainty of the astrometric calibration, propagated to any position on the image.

The covariance of the free platepar parameters (pointing and forward distortion) is estimated from the
residuals of the calibration stars stored in the platepar, and propagated to the sky direction of any pixel.
This is the error of the calibration itself, shared by all points measured with the same platepar. The
random error of each individual pick is not included.
"""

from __future__ import print_function, division, absolute_import

import argparse
import contextlib
import copy
import io
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from RMS.Formats.Platepar import Platepar, getPairedStarsSkyPositions


def freeParams(platepar, fixed_scale=False, only_pointing=False):
    """ Vector of the free parameters which map image coordinates to the sky. """

    params = [platepar.RA_d, platepar.dec_d, platepar.pos_angle_ref]

    if not fixed_scale:
        params.append(platepar.F_scale)

    if not only_pointing:
        params += list(platepar.x_poly_fwd[:platepar.poly_length])

        # Radial distortions store all their parameters in the X polynomial
        if platepar.distortion_type.startswith("poly"):
            params += list(platepar.y_poly_fwd[:platepar.poly_length])

    return np.array(params, dtype=np.float64)


def setFreeParams(platepar, params, fixed_scale=False, only_pointing=False):
    """ Copy of the platepar with the given free parameters (the inverse of freeParams). """

    pp = copy.copy(platepar)
    pp.RA_d, pp.dec_d, pp.pos_angle_ref = params[:3]
    i = 3

    if not fixed_scale:
        pp.F_scale = params[i]
        i += 1

    if not only_pointing:
        n = pp.poly_length

        pp.x_poly_fwd = np.array(pp.x_poly_fwd)
        pp.x_poly_fwd[:n] = params[i:i + n]

        if pp.distortion_type.startswith("poly"):
            pp.y_poly_fwd = np.array(pp.y_poly_fwd)
            pp.y_poly_fwd[:n] = params[i + n:i + 2*n]

    return pp


def skyOffsets(platepar, jd, x, y, ra_ref, dec_ref):
    """ Tangent plane offsets (arcsec) of the sky directions of pixels (x, y) from the reference directions,
        stacked as [offsets along RA*cos(dec), offsets along Dec].
    """

    ra, dec = getPairedStarsSkyPositions(x, y, jd, platepar)
    d_ra = (ra - ra_ref + 180)%360 - 180

    return 3600*np.concatenate([d_ra*np.cos(np.radians(dec_ref)), dec - dec_ref])


def numericalJacobian(func, params, target_shift=0.01):
    """ Central difference Jacobian of func. The parameters span many orders of magnitude (degrees for the
        pointing, polynomial coefficients down to ~1e-9), so every step is scaled to move the mapped
        positions by about target_shift arcsec.
    """

    f0 = func(params)
    jac = np.zeros((len(f0), len(params)))

    for j in range(len(params)):

        step = np.zeros_like(params)
        step[j] = 1e-8*max(1.0, abs(params[j]))
        step[j] *= target_shift/max(np.max(np.abs(func(params + step) - f0)), 1e-12)

        jac[:, j] = (func(params + step) - func(params - step))/(2*step[j])

    return jac


def plateparCovariance(platepar, fixed_scale=False, only_pointing=False):
    """ Covariance of the free platepar parameters, from the residuals of the calibration stars.

    Return:
        (params, cov, sigma, cond):
            - params: [ndarray] Free parameters (see freeParams).
            - cov: [ndarray] Their covariance matrix.
            - sigma: [float] Star residual per axis (arcsec), corrected for the number of fitted parameters.
            - cond: [float] Condition number of the normalized Jacobian (large = poorly constrained fit).
    """

    stars = np.array(platepar.star_list)
    jd, x, y, ra_cat, dec_cat = stars[0, 0], stars[:, 1], stars[:, 2], stars[:, 4], stars[:, 5]

    params = freeParams(platepar, fixed_scale, only_pointing)
    residuals = lambda p: skyOffsets(setFreeParams(platepar, p, fixed_scale, only_pointing), jd, x, y,
        ra_cat, dec_cat)

    res = residuals(params)
    jac = numericalJacobian(residuals, params)

    sigma = np.sqrt(np.sum(res**2)/(len(res) - len(params)))

    # cov = sigma^2*(J^T J)^-1, inverted through the SVD of the column-normalized Jacobian for stability
    norm = np.linalg.norm(jac, axis=0)
    u, s, vt = np.linalg.svd(jac/norm, full_matrices=False)
    cov = sigma**2*(vt.T/s**2) @ vt/np.outer(norm, norm)

    # At the least squares optimum the residuals are orthogonal to the Jacobian. If not, the platepar is not
    #   the best fit of its own model, and its real errors are larger than cov predicts
    reducible = np.sum((u.T @ res)**2)/np.sum(res**2)
    if reducible > 0.01:
        print("WARNING: the platepar is not at the optimum, refitting could remove {:.0%} of the squared "
            "residuals. The uncertainty is underestimated.".format(reducible))

    return params, cov, sigma, s[0]/s[-1]


def pointCovariance(platepar, params, cov, jd, x, y, fixed_scale=False, only_pointing=False):
    """ Propagate the parameter covariance to the sky directions of pixels (x, y).

    Return:
        (var_ra, var_dec, cov_ra_dec): [ndarrays] Covariance terms of every point (arcsec^2), RA along
            RA*cos(dec).
    """

    ra0, dec0 = getPairedStarsSkyPositions(x, y, jd, platepar)
    offsets = lambda p: skyOffsets(setFreeParams(platepar, p, fixed_scale, only_pointing), jd, x, y,
        ra0, dec0)

    jac = numericalJacobian(offsets, params)
    j_ra, j_dec = jac[:len(x)], jac[len(x):]

    var_ra = np.einsum('ij,jk,ik->i', j_ra, cov, j_ra)
    var_dec = np.einsum('ij,jk,ik->i', j_dec, cov, j_dec)
    cov_ra_dec = np.einsum('ij,jk,ik->i', j_ra, cov, j_dec)

    return var_ra, var_dec, cov_ra_dec


def monteCarloCovariance(platepar, sigma, x, y, n_runs, fixed_scale=False, only_pointing=False):
    """ Empirical check of pointCovariance: refit the platepar n_runs times to synthetic stars (the sky
        positions the platepar maps exactly, plus Gaussian noise of sigma arcsec per axis) and measure the
        mean squared error of the mapped pixels (x, y).
    """

    stars = np.array(platepar.star_list)
    jd, img_stars = stars[0, 0], stars[:, 1:4]

    ra_true, dec_true = getPairedStarsSkyPositions(stars[:, 1], stars[:, 2], jd, platepar)
    ra0, dec0 = getPairedStarsSkyPositions(x, y, jd, platepar)

    errors = []
    for _ in range(n_runs):

        noise_ra, noise_dec = np.random.normal(0, sigma/3600, (2, len(stars)))
        catalog = np.c_[ra_true + noise_ra/np.cos(np.radians(dec_true)), dec_true + noise_dec, stars[:, 6]]

        pp = copy.deepcopy(platepar)
        with contextlib.redirect_stdout(io.StringIO()):
            pp.fitAstrometry(jd, img_stars, catalog, fixed_scale=fixed_scale, fit_only_pointing=only_pointing)

        errors.append(skyOffsets(pp, jd, x, y, ra0, dec0))

    errors = np.array(errors)
    e_ra, e_dec = errors[:, :len(x)], errors[:, len(x):]

    return np.mean(e_ra**2, axis=0), np.mean(e_dec**2, axis=0), np.mean(e_ra*e_dec, axis=0)


def semiMajorAxis(var_ra, var_dec, cov_ra_dec):
    """ 1-sigma semi-major axis of the error ellipse. """

    return np.sqrt((var_ra + var_dec)/2 + np.sqrt(((var_ra - var_dec)/2)**2 + cov_ra_dec**2))



if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser(description="Map the uncertainty of the astrometric calibration over "
        "the image.")

    arg_parser.add_argument('platepar', type=str, help="Path to the platepar file.")

    arg_parser.add_argument('--fixedscale', action="store_true", help="The scale was kept fixed in the fit.")

    arg_parser.add_argument('--pointing', action="store_true", help="Only the pointing was fitted.")

    arg_parser.add_argument('--mc', type=int, default=0, metavar='N', help="Compare with N Monte Carlo refits.")

    cml_args = arg_parser.parse_args()


    pp = Platepar()
    pp.read(cml_args.platepar)

    if not pp.star_list:
        print("The platepar has no calibration stars stored.")
        sys.exit()

    stars = np.array(pp.star_list)
    jd = stars[0, 0]
    opts = dict(fixed_scale=cml_args.fixedscale, only_pointing=cml_args.pointing)

    params, cov, sigma, cond = plateparCovariance(pp, **opts)

    print("Stars: {:d}, free parameters: {:d} ({:s})".format(len(stars), len(params), pp.distortion_type))
    print("Star residual per axis: {:.1f} arcsec, condition number: {:.1e}".format(sigma, cond))


    # Grid of image positions with square cells
    nx = 64
    ny = int(round(nx*pp.Y_res/pp.X_res))
    grid_x, grid_y = np.meshgrid((np.arange(nx) + 0.5)*pp.X_res/nx, (np.arange(ny) + 0.5)*pp.Y_res/ny)
    x, y = grid_x.ravel(), grid_y.ravel()

    maps = [("Linear propagation", semiMajorAxis(*pointCovariance(pp, params, cov, jd, x, y, **opts)))]

    if cml_args.mc:
        maps.append(("Monte Carlo ({:d} refits)".format(cml_args.mc),
            semiMajorAxis(*monteCarloCovariance(pp, sigma, x, y, cml_args.mc, **opts))))

    for title, sig in maps:
        print("{:s}: 1-sigma semi-major axis min {:.1f}, median {:.1f}, max {:.1f} arcsec".format(title,
            np.min(sig), np.median(sig), np.max(sig)))

    if cml_args.mc:
        ratio = maps[1][1]/maps[0][1]
        print("Monte Carlo / linear: median {:.2f}, 5-95% {:.2f}-{:.2f}".format(np.median(ratio),
            *np.percentile(ratio, [5, 95])))


    # Plot the maps on a common logarithmic color scale, with the calibration stars
    fig, axes = plt.subplots(1, len(maps), figsize=(6*len(maps), 4), squeeze=False)
    norm = LogNorm(min(np.min(sig) for _, sig in maps), max(np.max(sig) for _, sig in maps))

    for ax, (title, sig) in zip(axes[0], maps):
        im = ax.imshow(sig.reshape(ny, nx), extent=[0, pp.X_res, pp.Y_res, 0], norm=norm)
        ax.scatter(stars[:, 1], stars[:, 2], s=6, c='w', edgecolors='k', linewidths=0.5)
        ax.set_title(title)

    fig.colorbar(im, ax=axes[0].tolist(), label="Calibration 1$\\sigma$ semi-major axis (arcsec)")

    plot_path = os.path.splitext(cml_args.platepar)[0] + "_uncertainty.png"
    plt.savefig(plot_path, dpi=150)
    print("Saved:", plot_path)

    plt.show()
