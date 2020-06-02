
path_to_dir=../preprint
source_folder=$path_to_dir/results/model_promoters_and_random_combined/

target_folder=... #accessible through internet




#bedGRaph to bed for mymethods


export PATH=$path_to_dir/softwares/utils/:$PATH


for cell_line in K562 GM12878
do
    echo $cell_line
    
    for distance_measure in ML Bayes_estimated_priors
    do
        echo $distance_measure
        predictions_path=$source_folder
        results_path=$predictions_path$cell_line"/"$distance_measure"/"


        
        mkdir $results_path"GenomeBrowser/"
	
        cd $results_path


 	      array=(*enhancers_PREPRINT*0.5.bedGraph)
	      for arr in "${array[@]}"
	      do
            awk '{OFS=" "; print $1, $2, $3}' $arr >  $results_path"GenomeBrowser/"${arr%"."*}".bed"
	          sort -k1,1 -k2,2n $results_path"GenomeBrowser/"${arr%"."*}".bed" > $results_path"GenomeBrowser/"${arr%"."*}".sorted.bed"
            bedToBigBed $results_path"GenomeBrowser/"${arr%"."*}".sorted.bed"  $path_to_dir/softwares/utils/hg19.chrom.sizes $results_path"GenomeBrowser/"${arr%"."*}".bb"
            rm $results_path"GenomeBrowser/"${arr%"."*}".sorted.bed"
 	      done


        array=(*prediction_scores*)
        sort -k1,1 -k2,2n $array > $results_path"GenomeBrowser/prediction_scores_PREPRINT.sorted.bedGraph"
        bedGraphToBigWig $results_path"GenomeBrowser/prediction_scores_PREPRINT.sorted.bedGraph" $path_to_dir/softwares/utils/hg19.chrom.sizes $results_path"GenomeBrowser/prediction_scores_PREPRINT.bw" 
        rm $results_path"GenomeBrowser/prediction_scores_PREPRINT.sorted.bedGraph"
    done
done

#check!


#RFECS

source_folder=$path_to_dir/results/


for cell_line in K562 GM12878
do
    results_path=$source_folder"RFECS_combined/"$cell_line"/bedfiles/"
    mkdir $results_path"GenomeBrowser/"


    for predictions_TSS in all without_TSS
    do

        arr="enhancers_"$predictions_TSS"_RFECS_threshold_05.bedGraph" 
	      awk '{OFS=" "; print $1, $2, $3}' $results_path$arr >  $results_path"GenomeBrowser/"${arr%"."*}".bed"
	      sort -k1,1 -k2,2n $results_path"GenomeBrowser/"${arr%"."*}".bed" > $results_path"GenomeBrowser/"${arr%"."*}".sorted.bed"
        bedToBigBed $results_path"GenomeBrowser/"${arr%"."*}".sorted.bed" $path_to_dir/softwares/utils/hg19.chrom.sizes $results_path"GenomeBrowser/"${arr%"."*}".bb"
        rm $results_path"GenomeBrowser/"${arr%"."*}".sorted.bed"
         
    done
  
    #all prediction scores for RFECS, blacklists removed from these
 
    array="all_prediction_scores_RFECS_threshold_05.bedGraph" 
    sort -k1,1 -k2,2n $results_path$array > $results_path"GenomeBrowser/"${array%"."*}".sorted.bedGraph"
    bedGraphToBigWig $results_path"GenomeBrowser/"${array%"."*}".sorted.bedGraph" $path_to_dir/softwares/utils/hg19.chrom.sizes $results_path"GenomeBrowser/"${array%"."*}".bw" 
    rm $results_path"GenomeBrowser/"${array%"."*}".sorted.bedGraph"    

done


#check!

##############Create direcs####################

#for random_str in pure_random random_with_signal
#do
random_str="combined"
mkdir $target_folder""$random_str
mkdir $target_folder""$random_str"/PREPRINT"
mkdir $target_folder""$random_str"/RFECS"
mkdir $target_folder""$random_str"/trackHub_GM12878"
mkdir $target_folder""$random_str"/trackHub_K562"
for cell_line in K562 GM12878
do

    mkdir $target_folder""$random_str"/PREPRINT/"$cell_line
    mkdir $target_folder""$random_str"/RFECS/"$cell_line

done
done

########Copy GenomeBrowser files to target!#############################################

#mypredictions

predictions_path=$path_to_dir/results/model_promoters_and_random_combined/



threshold=0.5
type=maxscore


random_str="combined"
for cell_line in K562 GM12878
do

    for distance_measure in ML Bayes_estimated_priors
    do

       
        results_path=$predictions_path$cell_line"/"$distance_measure"/GenomeBrowser/"

        
        
        
        for predictions_TSS in all without_TSS
       	do
	    
      	    scp $results_path"enhancers_PREPRINT_"$predictions_TSS"_"$type"_"$threshold".bb" $target_folder$random_str"/PREPRINT/"$cell_line"/"$distance_measure"_enhancers_"$predictions_TSS"_"$type"_"$threshold".bb"

	      done

        filename=$results_path"prediction_scores_PREPRINT.bw" 
        
	      scp $filename $target_folder$random_str"/PREPRINT/"$cell_line"/"$distance_measure"_prediction_scores_PREPRINT.bw"

    
    done
done
        
done
#check!

source_folder=$path_to_dir/results/




for cell_line in K562 GM12878
do

    
    results_path=$source_folder"/RFECS_combined/"$cell_line"/bedfiles/GenomeBrowser/"
    for predictions_TSS in all without_TSS
    do

        arr="enhancers_"$predictions_TSS"_RFECS_threshold_05.bed" 
	      scp $results_path""${arr%"."*}".bb" $target_folder$random_str"/RFECS/"$cell_line"/enhancers_"$predictions_TSS".bb"

         
    done
  
    #all prediction scores for RFECS, blacklists removed from these
 
    array="all_prediction_scores_RFECS_threshold_05" 
    scp $results_path""${array%"."*}".bw" $target_folder$random_str"/RFECS/"$cell_line"/all_prediction_scores.bw"


done

done

#check
