#!/usr/bin/awk -f
# A script that I got from the SILVA support to remove duplicate sequences from the SILVA SSU NR99 REF DB

{
	if ($0 ~ expression) {
   skip = !printMatch;
	 } else if (/^>/) {
	     skip = printMatch;
	 }

	 if (!skip) {
	     print;
	  }
}
