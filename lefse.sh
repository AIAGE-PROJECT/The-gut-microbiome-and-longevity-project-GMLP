#!/bin/bash
#Linear discriminant analysis Effect Size(lefse)
conda activate lefse
wd=/aiage/psf/project/guangxi_age1575/lefse/Figure4
lefse_dir=/aiage/software/lefse
cd $wd

#Fig4a-b and Supplementary Fig.5b
for i in level1-6_66-85vs100-117 level1-6_20-44vs100-117;
do python $lefse_dir/format_input.py $wd/data/${i}.txt $wd/${i}.in -c 2 -u 1 -o 100000;
   python $lefse_dir/run_lefse.py $wd/${i}.in $wd/${i}.res -l 2;
   python $lefse_dir/plot_res.py $wd/${i}.res $wd/${i}-lefse.pdf --format pdf --dpi 600 --feature_font_size 12 --class_legend_font_size 12 --left_space 0.2 --right_space 0.2 --width 10;
   python $lefse_dir/plot_cladogram.py $wd/${i}.res $wd/${i}-cladogram.pdf --format pdf  --dpi 600 --right_space_prop 0.3 --label_font_size 8 --class_legend_font_size 10;
done

#Fig4d and Supplementary Fig.5a
for i in level7_20-44_100-117vs66-85-subclass level1-6_20-44_100-117vs66-85-subclass;
do python $lefse_dir/format_input.py $wd/data/${i}.txt $wd/${i}.in -c 2 -s 3 -u 1 -o 100000;
   python $lefse_dir/run_lefse.py $wd/${i}.in $wd/${i}.res -l 2;
   python $lefse_dir/plot_res.py $wd/${i}.res $wd/${i}-lefse.pdf --format pdf --dpi 600 --feature_font_size 12 --class_legend_font_size 12 --left_space 0.2 --right_space 0.2 --width 12;
   python $lefse_dir/plot_cladogram.py $wd/${i}.res $wd/${i}-cladogram.pdf --format pdf  --dpi 600 --right_space_prop 0.3 --label_font_size 8 --class_legend_font_size 10;
done 
conda deactivate