/ *ERROR_SPRINTF.*/ {
if($0 == "#define ERROR_SPRINTF sprintf")
      nextfile
count = split($0,s,/[{]/)
i = 2;
while(i <= count) {
	out = out s[i]
	i = i+1
	}
enum = substr(out,match(out,"[0-9]{3}"),3)
if(para_match($0)!=0){
	getline
	i = match($0,"[a-zA-Z0-9]")
	out = out substr($0,i,length - i)
	out = substr(out,1,length(out)-1)
}
else
	out = substr(out,1,length(out)-2)
print "<li> " enum " \\ref " FILENAME "::" FNR " \"{" out > "./doc/text/background_errors_tmp.doc"
out = ""
}

function para_match(s) {
	 return (split(s,a,"(") - split(s,b,")"))
}
