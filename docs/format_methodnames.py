# format_methodnames.py
"""Do a find and replace operation on all code files in the repository.

Examples
--------

# Replace "Time_Instance_Objective" with "g"
$ python3 format_methodnames.py Time_Instance_Objective g
"""

import os

EXTENSIONS = (
    ".m",
    ".md",
)


def find_and_replace(filename: str, old: str, new: str) -> bool:
    """Do a find-and-replace in the given file.

    Parameters
    ----------
    filename : str
        Markdown file to read and format.
    old : str
        String to be replaced.
    new : str
        String to replace ``old`` with.

    Returns
    -------
    bool : True if the file changed, False otherwise.
    """
    # Read the file.
    with open(filename, "r") as infile:
        contents = infile.read()

    if old not in contents:
        return False

    # Write the file.
    with open(filename, "w") as outfile:
        outfile.write(contents.replace(old, new))

    return True


def format_all_files(
    directory: str,
    old: str,
    new: str,
    verbose: bool = False,
) -> bool:
    """Format all MATLAB code blocks in every Markdown (.md) file in the
    directory (descending through the directory tree recursively).

    Parameters
    ----------
    directory : str
        Directory to search for Markdown files to correct.
    old : str
        String to be replaced.
    new : str
        String to replace ``old`` with.
    verbose : bool
        If True, print before/after for each formatted block
        and report which files were changed.

    Returns
    -------
    bool : True if any files changed, False otherwise.
    """
    anyupdates = False
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if any(filename.endswith(ext) for ext in EXTENSIONS):
                thefile = os.path.join(dirpath, filename)
                updated = find_and_replace(thefile, old, new)
                if updated and verbose:
                    print(f"Updated {thefile}")
    return anyupdates


# =============================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    fname = os.path.basename(__file__)
    parser.usage = f""" python3 {fname} --help
        python3 {fname} OLD NEW [--quiet]
    """

    parser.add_argument(
        "old",
        type=str,
        help="string to replace",
    )
    parser.add_argument(
        "new",
        type=str,
        help="string to replace `old` with",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="do not print the list of altered files",
    )

    args = parser.parse_args()
    format_all_files("..", args.old, args.new, not args.quiet)
