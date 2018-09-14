#!/bin/bash

function fixname {
    index=`printf %04d $1`
    old=$2$number.root
    new=$3$index.root
    if [ -f $4/$new ]; then
        echo file $new already exists
        return 
    fi
    if [ -f $4/$old ]; then
        mv $4/$old $4/$new
        echo changed $old into $new in dir $4
    else 
        echo file $old not found
        return
    fi
}

for number in {1..500}; do
    fixname $number SimAnalyzed_Beta_ SimAnalyzed-beta-vector- $UCNA_MONTE_CARLO_DIR
    fixname $number SimAnalyzed_Beta_fierz_ SimAnalyzed-beta-fierz- $UCNA_MONTE_CARLO_DIR
    fixname $number SimAnalyzed_Beta_Fierz_ SimAnalyzed-beta-fierz- $UCNA_MONTE_CARLO_DIR
    fixname $number SimAnalyzed_Beta_twittled- SimAnalyzed-beta-vector-twittled- $UCNA_SYSTEMATICS_DIR
    fixname $number SimAnalyzed_Beta_fierz_twittled- SimAnalyzed-beta-fierz-twittled- $UCNA_SYSTEMATICS_DIR
    fixname $number SimAnalyzed_Beta_Fierz_twittled- SimAnalyzed-beta-fierz-twittled- $UCNA_SYSTEMATICS_DIR
done

