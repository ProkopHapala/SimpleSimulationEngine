import sys
import re

comment_pattern = re.compile(r'//.*?$|/\*.*?\*/', re.MULTILINE | re.DOTALL)

function_pattern = r'''
    # 1. Return type: alphanumeric with possible *
    ([a-zA-Z0-9_*]+)\s+
    # 2. Function name with mandatory space separation
    ([a-zA-Z0-9_]+)\s*
    # 3. Opening bracket with optional whitespace
    \(\s*
    # 4. Arguments - anything except brackets
    ([^{}()]*?)
    # 5. Closing bracket
    \)\s*
    # 6. Opening curly brace
    {
'''
function_pattern = re.compile( function_pattern, re.VERBOSE | re.MULTILINE | re.DOTALL )

def print_function_header( rec ):
    #f"{return_type} {func_name}({args})"
    return_type, func_name, args = rec
    print( f"{return_type} {func_name}({args})" )

def find_function_headers(filepath, bNoComments=True):
    with open(filepath, 'r') as f:
        content = f.read()
    if bNoComments:
            content = comment_pattern.sub('', content)
    matches = function_pattern.finditer(content)
    headers = []
    for match in matches:
        return_type, func_name, args = match.groups()
        header = ( return_type, func_name, args )
        headers.append(header)
    return headers

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python function_header_extractor.py <cpp_file>")
        exit()
        
    headers = find_function_headers(sys.argv[1])
    for h in headers:
        #print("\n" )
        print_function_header( h )