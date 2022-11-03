a=$(sed -n '3p' infile.in | cut -c2-6)
b=$(sed -n '4p' infile.in | cut -c2-6)
cat p.gnu | sed -e "s/AAA/$a/g" -e "s/BBB/$b/g"
