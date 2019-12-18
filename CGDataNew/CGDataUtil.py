import re

def col_fix( name ):
    newName=""
    regex = re.compile('[a-z,A-Z,0-9,_, ]')
    for i in range (0,len(name)):
        if not regex.match(name[i]):
            if i==0:
                continue
            if newName[-1]!="_":
                newName=newName+"_"
        else:
            newName =newName+name[i]
    if newName[-1]=="_":
        newName=newName[:-1]
    while (len(newName) > 64):
        out = re.sub( r'[aeiou]([^aioeu]*)$', r'\1', newName)
        if newName != out:
            newName =out
        else:
            newName = newName[:-1]
    return newName

def colStringLength (states):
    totalL =0
    collect =[]
    for state in states:
        if state == None:
            continue
        if state not in collect:
            collect.append(state)
            totalL = totalL +len(str(state))
    return totalL
    
def trackName_fix( name ):
    if trackName_good(name):
        return name
    if len(name) + len("genomic_") + len("_downsample") >=64:
        return False
    out=""
    regex = re.compile('[a-z,A-Z,0-9,_]')
    for i in range (0,len(name)):
        if not regex.match(name[i]):
            out= out+"_"
        else:
            out = out+name[i]
    return out

def trackName_good( name ):
    regex = re.compile('[a-z,A-Z,0-9,_]')
    for i in range (0,len(name)):
        if not regex.match(name[i]):
            return False
    if len(name) - len("genomic_") -len("_downsample")>64:
        return False
    return True


def sort_nicely( l ):
    """ Sort the given list in the way that humans expect.All items in l are expected to be string type
        e.g. x1, x2,x10
    """ 
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key ) 
    return
