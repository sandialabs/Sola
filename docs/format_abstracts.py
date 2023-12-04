# format_abstracts.py
"""Enable sphinxcontrib.matlab to document abstract MATLAB methods."""

import os
import re

# Path to MATLAB source code.
matlab_src_dir = os.path.abspath(os.path.join("..", "src"))

# Regular expressions.
_signature = r"^( *?)((?:\[?(?:\w|(?:, ))+?\]? = )?\w+?\(.*?\))"
_method_with_doc = re.compile(
    _signature + r"(\n\s+%.+?)\n\n",
    flags=(re.DOTALL | re.MULTILINE),
)
_method_without_doc = re.compile(
    _signature + r"\n\n",
    flags=re.MULTILINE,
)
_complete_method = re.compile(
    r"^ *?function.*?end",
    flags=(re.DOTALL | re.MULTILINE),
)

_hasdoc = (
    r"\1function \2\n"
    r"\1% *Abstract method.*\3\n"
    r"\1    error('not implemented');\n\1end\n\n"
)
_nodoc = r"\1function \2\n\1end\n\n"


def format_code_block(codeblock: str, verbose: bool = False):
    """Format a single code block.

    Parameters
    ----------
    codeblock : str
        Code block to format.
    verbose : bool
        If True, print information on what substitutions are made.

    Returns
    -------
    str
        Formatted code block.
    """
    # Remove "Abstract" from the "methods (...)" delimiter.
    matches = re.findall(
        r"methods\s*?\((.*?Abstract.*?)\)",
        codeblock,
        flags=re.DOTALL,
    )
    if not matches:
        raise RuntimeError("no 'methods (Abstract)' clause not found")
    if len(matches) > 1:
        raise ValueError("multiple 'methods (Abstract)' clauses found")

    preamble_args = [s.strip() for s in matches[0].split(", ")]
    preamble_args.remove("Abstract")
    preamble = ", ".join(preamble_args)
    if len(preamble) > 0:
        preamble = f" ({preamble})"
    codeblock = re.sub(
        r"methods\s*?\(.*?Abstract.*?\)",
        f"methods{preamble}",
        codeblock,
        flags=re.DOTALL,
    )

    # Format methods without a docstring.
    if verbose:
        for match in _method_without_doc.findall(codeblock):
            print(f"_method_without_doc: Modifying {match[1]}")
    codeblock = _method_without_doc.sub(_nodoc, codeblock)

    # Format methods with a docstring.
    if verbose:
        for match in _method_with_doc.findall(codeblock):
            print(f"_method_with_doc: Modifying {match[1]}")
            # for i, m in enumerate(match):
            #     print(f"match[{i}]: '{m}'")
            #     print()
    codeblock = _method_with_doc.sub(_hasdoc, codeblock)

    # Fix docstring indentation.
    for match in _complete_method.findall(codeblock):
        lines = match.split("\n")
        indent = " " * (len(re.findall(r"^(\s+?)function", lines[0])[0]) + 4)
        for i in range(1, len(lines)):
            if not (text := lines[i].strip()).startswith("%"):
                break
            lines[i] = indent + text
        codeblock = codeblock.replace(match, "\n".join(lines))

    return codeblock


def get_abstract_methods_blocks(code: str) -> list:
    """Extract blocks from the code starting with ``methods (Abstract)`` and
    ending with ``end``.

    Parameters
    ----------
    code : str
        Source code to examine (as a string).

    Returns
    -------
    blocks : list
        Abstract methods blocks.
    """
    return re.findall(
        r"^\s*?methods\s*?\(.*?Abstract.*?\).+?end",
        code,
        flags=(re.DOTALL | re.MULTILINE),
    )


def format_file(filename: str) -> None:
    """Format the abstract code blocks of the specified file."""
    with open(filename, "r") as infile:
        code = infile.read()

    for block in get_abstract_methods_blocks(code):
        newblock = format_code_block(block)
        if newblock != block:
            print(f"Making a change to {os.path.basename(filename)}")
            code = code.replace(block, newblock, 1)

    with open(filename, "w") as outfile:
        outfile.write(code)


def main():
    """Walk though the source directory and format all of the MATLAB files."""
    for dirpath, _, filenames in os.walk(matlab_src_dir):
        for filename in filenames:
            if filename.endswith(".m"):
                format_file(os.path.join(dirpath, filename))


# =============================================================================

_test_in = r"""
methods    ( Abstract, Access = public, Static)

    abstractmethod1(arg1)
    % Docstring for abstractmethod1().
    %
    % Parameters
    % ----------
    % arg1
    %   First argument.

    [out1, out2] = abstractmethod2()
    % Docstring for abstractmethod2().
    %
    % Returns
    % -------
    % out1
    %   First return.
    % out2
    %   Second return.

    out = abstractmethod3(arg1, arg2)
    % Short docstring.

    abstractmethodwithnodocstring()

    thing = doclessmethod(arg1)

end
"""

_test_out = r"""
methods (Access = public, Static)

    function abstractmethod1(arg1)
        % *Abstract method.*
        % Docstring for abstractmethod1().
        %
        % Parameters
        % ----------
        % arg1
        %   First argument.
        error('not implemented');
    end

    function [out1, out2] = abstractmethod2()
        % *Abstract method.*
        % Docstring for abstractmethod2().
        %
        % Returns
        % -------
        % out1
        %   First return.
        % out2
        %   Second return.
        error('not implemented');
    end

    function out = abstractmethod3(arg1, arg2)
        % *Abstract method.*
        % Short docstring.
        error('not implemented');
    end

    function abstractmethodwithnodocstring()
    end

    function thing = doclessmethod(arg1)
    end

end
"""

_test_in2 = r"""
classdef My_Class < handle

    methods (Abstract, Silly)

        method1()

        method2(arg1)

        out = method3(arg1, arg2)
            % This one has a docstring.

        [out1, out2] = method3(arg1, arg2)
            % as does this one.
    end

    methods (Abstract)

        method5()

    end

    methods

        function method5()
        end

    end

end
"""

_test_out2_0 = r"""
    methods (Abstract, Silly)

        method1()

        method2(arg1)

        out = method3(arg1, arg2)
            % This one has a docstring.

        [out1, out2] = method3(arg1, arg2)
            % as does this one.
    end"""

_test_out2_1 = r"""
    methods (Abstract)

        method5()

    end"""

__doc__ = f"""
Enable sphinxcontrib.matlab to document abstract MATLAB methods
by temporarily rewriting them as non-abstract methods.

This is step 2 in the following sequence.

1. Switch to a new, temporary git branch.
2. Comb through the MATLAB source code and replace abstract methods
   with non-abstract methods.
3. Build the documentation.
4. Switch back to the old branch.
5. Delete the temporary branch.

Example
-------

Consider the following MATLAB code block.

```matlab
{_test_in}
```

This is replaced with the following.

```matlab
{_test_out}
```
"""


def _tests():
    """Unit tests."""
    # Test get_abstract_methods_block().
    outputs = get_abstract_methods_blocks(_test_in2)
    if len(outputs) != 2:
        raise AssertionError("get_abstract_methods_blocks() failed!")
    for i, out in enumerate([_test_out2_0, _test_out2_1]):
        if outputs[i] != out:
            print(
                f"Expected\n--------\n'{out}'\n"
                f"\nReceived\n--------\n'{outputs[i]}'"
            )
            raise AssertionError("get_abstract_methods_blocks() failed!")

    # Test format_code_block()
    if (output := format_code_block(_test_in, verbose=False)) != _test_out:
        format_code_block(_test_in, verbose=True)
        print(f"Input\n-----\n{_test_in}\n\nOutput\n------\n{output}")
        raise AssertionError("format_code_block() failed!")

    print("All tests passed!")


# =============================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    fname = os.path.basename(__file__)
    parser.usage = f""" python3 {fname} --help
        python3 {fname} --test
        python3 {fname}
    """

    parser.add_argument("--test", action="store_true", help="run unit tests")

    args = parser.parse_args()
    if args.test:
        _tests()
    else:
        main()
