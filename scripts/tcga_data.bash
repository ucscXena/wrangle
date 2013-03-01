export PATH=/inside/home/common/bin:$PATH; python -V 2>pythonVersion; which python >>pytonVersion;

cd /inside/home/jzhu/cgDataJing/TCGAscripts;
python run.py configCLIN logCLIN > changeCLIN;
grep "Feature Not in dictionary" logCLIN |sort |uniq > missing;
grep "feature name is modified" changeCLIN  |sort |uniq > change;
python featureSort.py;
mv TCGAfeature.py TCGAfeature.py_BK;
mv missingpython TCGAfeature.py;

python runBATCH.py configBATCH log;

cd /inside/home/jzhu/cgDataJing/scripts;
python runFlattenTCGA.py;

cd /inside/home/jzhu/cgDataUCSC2/scripts;
python compileBulkTCGA.py;
