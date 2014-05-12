/\\begin\{essyntax\}/ { inenv=1; print $0; next }
/\\end\{essyntax\}/ { inenv=0; print $0; next }
inenv==1 { print $0 }
