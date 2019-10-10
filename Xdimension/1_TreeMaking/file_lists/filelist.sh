#! /bin/bash   

#list1=`samweb list-definition-files prod_crossingmu_right_SCE_mcc9.1_reco2`
#list1+=`samweb list-definition-files prod_crossingmu_left_SCE_mcc9.1_reco2`
#list1+=`samweb list-definition-files prod_crossingmu_topleft_SCE_mcc9.1_reco2`
#list1+=`samweb list-definition-files prod_crossingmu_bottomleft_SCE_mcc9.1_reco2`
#list1+=`samweb list-definition-files prod_crossingmu_topright_SCE_mcc9.1_reco2`
#list1+=`samweb list-definition-files prod_crossingmu_bottomright_SCE_mcc9.1_reco2`

list1=`samweb list-files "defname: yatesla_intimeACPTmuon_uboone_overlay_DLup_run1_reco2"`

i=1
for file in $list1
do
    fullpath=`samweb locate-file $file`
    cutpath=`echo "${fullpath%(*}"/$file`
    echo "${cutpath##*:}"
    let i=i+1
done