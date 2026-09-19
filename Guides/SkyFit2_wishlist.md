# SkyFit2 wishlist: status

Review of all branches, commits, and the 700 PRs and issues of CroatianMeteorNetwork/RMS (September 2026, base
`prerelease` a81cbf8d).

| Request | Status | Notes |
|---|---|---|
| Frame slider / go to frame N in manual reduction | Not done | See below |
| Names of bright stars | Done (prerelease) | "Show Star Names" checkbox, spectral type, click opens SIMBAD (PR #762) |
| Frame number in the ECSV | Not done | See below |
| `.state` saved in the working directory | Not done, by design | See below |
| Clicking a point in the astrometry plot highlights it on the image | Done (prerelease) | `onAstrometryPlotPick`, commit a6f13a06 |
| Calibrate a video with a video/image from another date | Partial | A platepar from another folder can be loaded; no merging of stars from two dates |
| Astrometric error for every point | Not done | Prototype in this branch, see `AstrometryUncertainty.md` |

## Frame slider

No technical reason was found against it:

- The only slider (`image_navigation_slider`) navigates FF files in SkyFit mode, not frames.
- Every video frame is already read with a random seek (`cap.set(1, frame)`), so jumping to any frame costs the same
  as one arrow key, and memory does not change. `nextImg(n)` already accepts any jump.
- Keys: CTRL + arrows jump 10 frames, up/down 25.

To keep in mind: DFN and single image modes only allow moving next to the picks, so the slider must be disabled
there; while dragging, jump only on release (H.264 seeks decode from the previous keyframe).

## Frame number in the ECSV

The ECSV follows the common (GDEF) format, where time is the reference. Adding a column is safe: both the SkyFit2 and
the WMPL readers find columns by name. To decide: store the original frame index (what the user sees), not the
time-normalized/rolling-shutter-corrected number (already in `FTPdetectinfo_*_manual.txt`).

## Location of the .state file

The `.state` is meant to live next to the data: on load, all paths (data, input, platepar) are rebuilt from the
folder of the `.state`, so a folder can be moved between computers. Saving it elsewhere would break loading, unless
the absolute data path is also stored and used as a fallback.

The name is always `skyFitMR_latest.state`, so every save overwrites the previous one without a backup.
