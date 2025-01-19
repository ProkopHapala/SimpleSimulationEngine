import re
import sys

# Pre-compile patterns as before
COMMENT_PATTERN     = re.compile(r'//.*?$|/\*.*?\*/', re.MULTILINE | re.DOTALL)
FUNCTION_PATTERN    = re.compile(r'''([a-zA-Z0-9_*]+)\s+([a-zA-Z0-9_]+)\s*\(\s*([^{}()]*?)\)\s*{''', re.VERBOSE | re.MULTILINE | re.DOTALL)
SCOPE_START_PATTERN = re.compile(r'\b(class|struct|namespace)\s+([a-zA-Z0-9_]+)')

CPP_KEYWORDS = {
    'if', 'else', 'while', 'for', 'switch', 'case', 'return',
    'break', 'continue', 'do', 'sizeof', 'typedef'
}

def is_not_function(return_type, name):
    return (  return_type.strip() in CPP_KEYWORDS 
           or name       .strip() in CPP_KEYWORDS)

class Scope:
    def __init__(self, type_, name, start_pos):
        self.type = type_
        self.name = name
        self.start_pos = start_pos

def print_function_header( fh ):
    if f['scope']:
        print(f"{fh['return_type']} {fh['scope']}::{fh['name']}({fh['args']})")
    else:
        print(f"{fh['return_type']} {fh['name']}({fh['args']})")

def analyze_scopes_and_functions(content, bPrint=False ):
    scope_stack = []
    functions = []
    
    pos = 0
    while pos < len(content):
        if content[pos] == '{':
            # Look back to determine scope type
            look_back = content[max(0, pos-100):pos]
            scope_match = SCOPE_START_PATTERN.search(look_back)
            
            if scope_match:
                scope_type, scope_name = scope_match.groups()
                scope = Scope(scope_type, scope_name, pos)
                if bPrint: print( (" "*4*len(scope_stack))+scope.type+" "+scope.name+"{" ) 
            else:
                scope = Scope('block', '', pos)
            
            scope_stack.append(scope)  
        elif content[pos] == '}':
            if scope_stack:
                if bPrint and scope_stack[-1].type != 'block':
                    print( (" "*4*(len(scope_stack)-1))+"}")
                scope_stack.pop()

                
        # Check for function matches at current position
        func_match = FUNCTION_PATTERN.match(content, pos)
        if func_match:
            return_type, func_name, args = func_match.groups()
            scope_path = '::'.join(s.name for s in scope_stack if s.type in ('class', 'struct', 'namespace') and s.name)
            #full_name = f"{scope_path}::{func_name}" if scope_path else func_name
            if not is_not_function(return_type, func_name):
                #print(f"{scope_path}::{func_name}({args})")
                if bPrint: print( (" "*4*len(scope_stack)) + f"{return_type} {func_name}({args})")
                functions.append({
                    'name': func_name,
                    'return_type': return_type,
                    'args': args,
                    'scope': scope_path
                })
            pos = func_match.end()-1
        else:
            pos += 1
        
        
    return functions

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python function_header_extractor.py <cpp_file>")
        exit()
    with open(sys.argv[1], 'r') as f: content = f.read()
    content   = COMMENT_PATTERN.sub('', content)
    functions = analyze_scopes_and_functions(content, True)
    print("\n ============ \n ")
    for f in functions:
        #print()
        #print(f"@ {f['scope']}")
        #print(f"{f['return_type']} {f['name']}({f['args']})")
        #print(f"{f['scope']}::{f['return_type']} {f['name']}({f['args']})")
        print_function_header( f )
