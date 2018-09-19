#!/bin/bash

function replaceinstring {
    if [[ "$#" -ne "3" ]] ; then 
        BenchPrintError "Error you must pass only three arguments to BenchReplaceInString ($# passed - $@ - ${FUNCNAME})"
    fi

    str="$1"
    while [[ `echo "$str" | grep "$2" | wc -l` -ne 0 ]] ; do
        str="${str/$2/$3}"
    done
    echo "$str"
}

function replaceinfile {
    if [[ "$VERBOSE" == "1" ]] ; then
        echo replaceinfile "$1" "$2" "$3"
    fi
    allfile=`cat "$1"`
    allfile=`replaceinstring "$allfile" "$2" "$3"`
    echo "$allfile" > "$1"
}

function insertbefore {
    if [[ "$VERBOSE" == "1" ]] ; then
        echo insertbefore "$1" "$2" "$3"
    fi

    posparse=`grep -n "$2" "$1" | grep -Eo '^[^:]+' `
    if [[ "$posparse" == "" ]] ; then
        echo "[ERROR] cannot find line $2 in file $1"
        exit 1
    fi
    output=` awk -v n=$posparse -v s="$3" 'NR == n {print s} {print}' "$1" `
    echo "$output" > "$1"
}

echo "======================================"
echo "Auto add domain tool."
echo "======================================"

if [[ ! -f "./add-domain.sh" ]] ; then
    echo "[ERROR] the script must be run its directory"
    exit 1
fi

metaStartParse="// AUTO-ADD-DOMAIN START-PARSE"
metaEndParse="// AUTO-ADD-DOMAIN END-PARSE"
metaParseCode="  else if (type == domains::@CLASSNAME@::Type()) {
    return parseGeneric<domains::@CLASSNAME@>(domainTypename, node, beadTypes, connections, id);
  }"

metaIncludeCode="#include \"domains/@FILENAME@\""
metaHeaderPos="// AUTO-ADD-DOMAIN HEADER"

metaRestoreCode="  else if (dom.type() == domains::@CLASSNAME@::Type()) {
    static_assert(false,
        \"Error, you must decide if the domain @STRINGDOMAIN@\"
        \"can use a generic method (as restoreGenericDom) to be restored.\"
        \"Generic methods only manage Domain's attributes.\");
    // Maybe restoreGenericDom(ncDomain, dom);
  }"
metaRestorePos="// AUTO-ADD-DOMAIN RESTORE-DOM"

metaSaveCode="  else if (dom.type() == domains::@CLASSNAME@::Type()) {
    static_assert(false,
        \"Error, you must decide if the domain @STRINGDOMAIN@\"
        \"can use a generic method (as saveGenericDom) to be saved.\"
        \"Generic methods only manage Domain's attributes.\");
    // Maybe saveGenericDom(group, dom);
  }"
metaSavePos="// AUTO-ADD-DOMAIN SAVE-DOM"

##########################################################
# Get domain class name
echo "Please, enter the domain class name."
echo "The name must start by a letter and use only letters or numeric caractere"
echo "Example: MyNewDomain"
echo -n ">> "
read type

if [[ "$type" == "q" ]] || [[ "$type" == "Q" ]] ; then
    exit 0
fi

if [[ "$type" == "" ]] ; then
    echo "[ERROR] domain name cannot be null"
    exit 1
fi

typeLower=` echo "$type" | tr '[:upper:]' '[:lower:]' `
typeUpper=` echo "$type" | tr '[:lower:]' '[:upper:]' `

##########################################################
echo "Please, enter the domain filename."
echo "Default is: $typeLower.h"
echo -n ">> "
read classFilename

if [[ "$classFilename" == "q" ]] || [[ "$classFilename" == "Q" ]] ; then
    exit 0
fi

if [[ "$classFilename" == "" ]] ; then
    classFilename="$typeLower.h"
fi

##########################################################
echo "Please, enter the domain config name (use as an id and in the config files)."
echo "Default is: $typeLower"
echo -n ">> "
read configName

if [[ "$configName" == "q" ]] || [[ "$configName" == "Q" ]] ; then
    exit 1
fi

if [[ "$configName" == "" ]] ; then
    configName="$typeLower"
fi

##########################################################
# Move

templateFile="./template-domain.h"
if [[ ! -f "$templateFile" ]] ; then
        echo "[ERROR] cannot find the template file $templateFile"
    exit 1
fi

domainFile="../../src/domains/$classFilename"

if [[ -f "$domainFile" ]] ; then
    echo "[ERROR] the destination file $domainFile already exist"
    exit 1
fi

if [[ "$VERBOSE" == "1" ]] ; then
    echo cp "$templateFile" "$domainFile"
fi

cp "$templateFile" "$domainFile"
if [[ "$?" != "0" ]] ; then
    echo "[ERROR] cp $templateFile $domainFile failled with error $?"
    exit 1
fi

##########################################################
# Replace

replaceinfile "$domainFile" "@UPPERCASEDOMAIN@" "$typeUpper"
replaceinfile "$domainFile" "@LOWERCASEDOMAIN@" "$typeLower"
replaceinfile "$domainFile" "@STRINGDOMAIN@" "$configName"
replaceinfile "$domainFile" "@CLASSNAME@" "$type"
replaceinfile "$domainFile" "@FILENAME@" "$classFilename"
replaceinfile "$domainFile" "@USER@" "$USER"
replaceinfile "$domainFile" "@DATE@" " `date` "

##########################################################

parserfile="../../src/io/cplx.cpp"

if [[ ! -f "$parserfile" ]] ; then
        echo "[ERROR] cannot find the source file $parserfile"
    exit 1
fi

insertbefore "$parserfile" "$metaEndParse" "$metaParseCode"
insertbefore "$parserfile" "$metaHeaderPos" "$metaIncludeCode"

replaceinfile "$parserfile" "@UPPERCASEDOMAIN@" "$typeUpper"
replaceinfile "$parserfile" "@LOWERCASEDOMAIN@" "$typeLower"
replaceinfile "$parserfile" "@STRINGDOMAIN@" "$configName"
replaceinfile "$parserfile" "@CLASSNAME@" "$type"
replaceinfile "$parserfile" "@FILENAME@" "$classFilename"

##########################################################

echo "The domain has been added."
echo "You need to use the keyword \"$configName\" in the configuration."
echo "You must edit your domain file: \"src/domains/$classFilename\""
echo "You may need to update \"src/io/netcdf.cpp\" to decide how the domain can be save and restore"

