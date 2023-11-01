# format_code_blocks.py
"""Format MATLAB code blocks within the documentation using mh_style."""

import os
import re
import sys
import subprocess


__TEMP_FILE = "__matlabcodeblock.m"


def format_matlab_code_block(code_block: str) -> str:
    """Format a single block of code.

    1) Write the code block to a temporary file.
    2) Format the file with the command-line tool `mh_style`.
    3) Read back the formatted file.
    4) Delete the intermediate file.

    Parameters
    ----------
    code_block : str
        MATLAB code to format (as a string).

    Returns
    -------
    str : Formatted code.
    """
    # Write the code block to a temporary file.
    with open(__TEMP_FILE, "w") as outfile:
        outfile.write(code_block)

    # Format the file with the command-line tool `mh_style`.
    subprocess.run(
        ["mh_style", "--brief", "--fix", __TEMP_FILE],
        stdout=subprocess.DEVNULL,
    )

    # Read back the formatted file.
    with open(__TEMP_FILE, "r") as infile:
        formatted_code = infile.read()

    # Delete the intermediate file.
    os.remove(__TEMP_FILE)

    return formatted_code


def format_markdown_file(filename: str, verbose: bool = False) -> bool:
    """Format all of the code blocks in the given Markdown file.

    Parameters
    ----------
    filename : str
        Markdown file to read and format.
    verbose : bool
        If True, print before/after for each formatted block.

    Returns
    -------
    bool : True if the file changed, False otherwise.
    """
    if not filename.endswith(".md"):
        raise ValueError("this method is for Markdown files only")

    # Read the file.
    with open(filename, "r") as infile:
        contents = infile.read()

    # Format each MATLAB code block.
    changed = False
    for match in re.findall(r"```matlab\n(.*?)```", contents, flags=re.DOTALL):
        if (formatted := format_matlab_code_block(match)) != match:
            if verbose:
                print(f"\nUnformatted\n-----------\n{match}")
                print(f"\nFormatted\n---------\n{formatted}")
            contents = contents.replace(match, f"{formatted}", 1)
            changed = True

    # Write the results back to the file.
    if changed:
        if verbose:
            print(f"Updating {filename}")
        with open(filename, "w") as outfile:
            outfile.write(contents)

    return changed


def format_all_files(directory: str, verbose: bool = False) -> bool:
    """Format all MATLAB code blocks in every Markdown (.md) file in the
    directory (descending through the directory tree recursively).

    Parameters
    ----------
    directory : str
        Directory to search for Markdown files to correct.
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
            if filename.endswith(".md"):
                anyupdates |= format_markdown_file(
                    os.path.join(dirpath, filename),
                    verbose=verbose,
                )
            # elif filename.endswith(".rst"):
            #     format_restructuredtext_file(  # TODO if needed
            #         os.path.join(dirpath, filename),
            #         verbose=verbose,
            #     )
    return anyupdates


if __name__ == "__main__":
    sys.exit(
        int(format_all_files(os.getcwd(), verbose=("--verbose" in sys.argv)))
    )
