/ *ERROR_SPRINTF.*/ {
  if($0 == "#define ERROR_SPRINTF sprintf")
    nextfile
  out = ""
  count = split($0,s,/{/)
  c = 2;
  while(c <= count){
    out = out s[c]
    c = c+1
  }
  errcode = substr(out,match(out,/[0-9]{3}/),3)
  if(para_match($0) != 0){
    getline
    c = match($0,/[a-zA-Z0-9]/)
    out = out substr($0,c,length-c+1)
    out = substr(out,1,length(out)-2)
  }
  else
    out = substr(out,1,length(out)-2)
  print "<li> " errcode " \\ref " FILENAME "::" fname " : \" {" out > "./doc/text/background_errors_tmp.doc"
}

/.*/ {
  if(FNR == 1){
    i = 0
    j = 0
    cflag = 0
  }
  nocomments()
  i = i + split(" " $0 " ",tmp,/{/)
  i = i - split(" " $0 " ",tmp,/}/)
  if(i == 1 && j == 0) {
    if(match($0,/[a-zA-Z0-9_]+\(.*\).*{/) == 0)
      $0 = buf $0
    x = match($0,/\(/)
    $0 = substr($0,1,x-1)
    x = split($0,s)
    fname = s[x]
    if(match(fname,/\*/) == 1)
      fname = substr(fname,2,length(fname)-1)
  }
  buf = $0
  j = i
}

function nocomments(){
  tmp0 = "";
  spos = match($0,/\/\*/)
  eposold = 0
  epos = match($0,/\*\//)
  while(spos > 0 || epos > 0){
    if((spos < epos && spos > 0) || epos == 0){
      if(cflag == 0){
        tmp0 = tmp0 substr($0,eposold+1,spos-eposold-1)
        cflag = 1
      }
    sub(/\/\*/,"XX",$0)
    spos = match($0,/\/\*/)
    }
    else{
      if(cflag == 1){
        cflag = 0
        eposold = epos+1
      }
      sub(/\*\//,"YY",$0)
      epos = match($0,/\*\//)
    }
  }
  if(cflag == 0)
    tmp0 = tmp0 substr($0,eposold+1,length-eposold)
  $0 = tmp0
}

function para_match(s){
	 return (split(s,a,"(") - split(s,b,")"))
}