"""Utility function to help find image data on Bandicoot (storage mnt) across multiple users and directories.
This function is adapted from Mike Lippincott and original code can be found here:
https://github.com/WayScience/NF1_3D_organoid_profiling_pipeline/blob/76e313d81ef0c4b60ed05228e31a15c95c0a1dfa/utils/src/image_analysis_3D/file_utils/notebook_init_utils.py#L77"""

import os
import pathlib

def bandicoot_check(
    bandicoot_mount_path: pathlib.Path, root_dir: pathlib.Path
) -> pathlib.Path:
    """
    This function determines if the external mount point for Bandicoot exists.

    Parameters
    ----------
    bandicoot_mount_path : pathlib.Path
        The path to the Bandicoot mount point.
    root_dir : pathlib.Path
        The root directory of the Git repository.

    Returns
    -------
    pathlib.Path
        The base directory for image data.
    """
    if bandicoot_mount_path.exists():
        # comment out depending on whose computer you are on
        # mike's computer
        image_base_dir = pathlib.Path(
            os.path.expanduser("~/mnt/bandicoot/")
        ).resolve(strict=True)
    else:
        image_base_dir = root_dir
    return image_base_dir
