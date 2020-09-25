#!/bin/bash

code=$1
anno_dir=$PWD/4_Annotation
features=$anno_dir/features/features_report.txt

ko2code=$anno_dir/ko2code_tmp
code2pathway=$anno_dir/code2pathway_tmp
pathway2img=$anno_dir/KEGG/Pathway_${code}

if [ $? -eq 0 ]; then
	echo "Start KO -> ${code}_ID -> ${code}_pathwayID -> pathway.png"
	
	cat $features | cut -f8 | sed 's/\.//g' | sed '/^\s*$/d' | sed 's/`/\n/g' | sed 1d | sed 's/KO://g' > ${anno_dir}/ko.list

	if [ ! -e ${anno_dir}/ko2code.list ]; then 
		if [ ! -d $ko2code ]; then
			mkdir $ko2code
		fi
		for ko in `cat ${anno_dir}/ko.list`; do
			if [ $? != 0 ]; then
echo "### WARNING: this step might be failed frequently due to curl.
If failed please try to re-run following script standalone more times!!!
			
% /export/EC1680U/perl/bin/Trinity/ko2pathway.sh $code"
				exit $?
			fi
			if [ ! -f ${ko2code}/${ko}.tmp ]; then
				echo "$ko ---> ${code}_ID..."
				curl -# http://rest.kegg.jp/link/${code}/$ko > ${ko2code}/${ko}.tmp
			fi
		done
		cat ${ko2code}/*.tmp > ${anno_dir}/ko2code.list
		rm -rf ${ko2code}
	fi

	if [ ! -e ${anno_dir}/code2pathway.list ]; then
		if [ ! -d $code2pathway ]; then
			mkdir $code2pathway
		fi
		for codeID in `cat ${anno_dir}/ko2code.list | cut -f2`; do
			if [ $? != 0 ]; then
echo "### WARNING: this step might be failed frequently due to curl.
If failed please try to re-run following script standalone more times!!!

% /export/EC1680U/perl/bin/Trinity/ko2pathway.sh $code"
				exit $?
			fi
			if [ ! -f ${code2pathway}/${codeID}.tmp ]; then
				echo "${codeID} ---> ${code}_pathwayID..."
				curl -# http://rest.kegg.jp/link/pathway/${codeID} > ${code2pathway}/${codeID}.tmp
			fi
		done
		cat ${code2pathway}/*.tmp > ${anno_dir}/code2pathway.list
		rm -rf ${code2pathway}
	fi

	if [ ! -e ${anno_dir}/${pathway2img}/pathway.ok ]; then
		if [ ! -d $pathway2img ]; then
    		mkdir $pathway2img
		fi
        for pathway in `cat ${anno_dir}/code2pathway.list | cut -f2 | awk -F':' '{print $2}'`; do
			if [ $? != 0 ]; then
echo "### WARNING: this step might be failed frequently due to curl.
If failed please try to re-run following script standalone more times!!!

% /export/EC1680U/perl/bin/Trinity/ko2pathway.sh $code"
				exit $?
			fi
			if [ ! -f ${pathway2img}/${pathway}.png ]; then
            	echo "${pathway} ---> pathway.png..."
            	curl -# http://rest.kegg.jp/get/${pathway}/image > ${pathway2img}/${pathway}.png
			fi
        done
		touch ${pathway2img}/pathway.ok
    fi
	echo "Complete!"
	exit 0
else 
	echo "There was a problem in data retrieval"
	exit 1
fi

