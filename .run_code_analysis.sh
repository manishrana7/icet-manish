#! /bin/bash
echo "Run pylint3: pylint3 icetdev/ | grep -s rated | awk '{print $7}' > .pylint_score "
pylint3 icetdev/ | grep -s rated | awk '{print $7}' > .pylint_score 
echo "Done."
echo ""
echo -n "pylint score: " 
cat .pylint_score

echo ""
