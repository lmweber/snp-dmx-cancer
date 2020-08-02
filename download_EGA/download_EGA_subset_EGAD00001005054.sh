# -------------------------------------------------------
# Script to download lung adenocarcinoma dataset from EGA
# -------------------------------------------------------

# downloading only intermediate and high TMB samples


# sample T08 (int TMB)
mkdir LUNG_T08
cd LUNG_T08
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482043
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482044
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482045
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482046
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482047
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482048
cd ..
# flatten directories
find LUNG_T08 -mindepth 2 -type f -exec mv -t LUNG_T08 -i '{}' +


# sample T09 (high TMB)
mkdir LUNG_T09
cd LUNG_T09
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482049
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482050
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482051
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482052
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482053
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482054
cd ..
# flatten directories
find LUNG_T09 -mindepth 2 -type f -exec mv -t LUNG_T09 -i '{}' +


# sample T20 (high TMB)
mkdir LUNG_T20
cd LUNG_T20
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482067
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482068
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482069
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482070
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482071
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482072
cd ..
# flatten directories
find LUNG_T20 -mindepth 2 -type f -exec mv -t LUNG_T20 -i '{}' +


# sample T25 (int TMB)
mkdir LUNG_T25
cd LUNG_T25
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482073
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482074
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482075
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482076
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482077
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482078
cd ..
# flatten directories
find LUNG_T25 -mindepth 2 -type f -exec mv -t LUNG_T25 -i '{}' +


# sample T28 (high TMB)
mkdir LUNG_T28
cd LUNG_T28
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482079
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482080
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482081
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482082
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482083
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482084
cd ..
# flatten directories
find LUNG_T28 -mindepth 2 -type f -exec mv -t LUNG_T28 -i '{}' +


# sample T31 (int TMB)
mkdir LUNG_T31
cd LUNG_T31
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482091
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482092
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482093
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482094
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482095
pyega3 -cf ~/EGA/credentials.json -c 50 fetch EGAF00002482096
cd ..
# flatten directories
find LUNG_T31 -mindepth 2 -type f -exec mv -t LUNG_T31 -i '{}' +

