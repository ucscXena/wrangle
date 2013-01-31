import string
import TCGAUtil

def longTitle():
    index= TCGAUtil.featureLongTitle.keys()
    index.sort()
    for key in index:
        print "\t\""+key+"\":\""+TCGAUtil.featureLongTitle[key]+"\","

def shortTitle():
    index= TCGAUtil.featureShortTitle.keys()
    index.sort()
    for key in index:
        print "\t\""+key+"\":\""+TCGAUtil.featureShortTitle[key]+"\","

longTitle()
print
#shortTitle()
