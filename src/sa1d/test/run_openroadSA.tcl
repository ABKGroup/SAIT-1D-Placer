set top_module ArtNet

set stamp1 [clock seconds]

set lef_dir "/home/memzfs_projects/ANG2.0/shyeon/pdk/asap7/lef"
read_lef ${lef_dir}/asap7_tech_1x_201209.lef
read_lef ${lef_dir}/asap7sc7p5t_27_L_1x_201211.lef
read_lef ${lef_dir}/asap7sc7p5t_27_R_1x_201211.lef
read_lef ${lef_dir}/asap7sc7p5t_27_SL_1x_201211.lef

set testcase $::env(TESTCASE)

if { $testcase == 1 } {
  read_def /home/memzfs_projects/ANG2.0/d3yoon/1Dtest/1stTestInvs/ArtNet_place.def
} elseif { $testcase == 2 } {
  read_def /home/memzfs_projects/ANG2.0/d3yoon/1Dtest/2ndTestInvs/ArtNet_place.def
} elseif { $testcase == 3 } {
  read_def /home/memzfs_projects/ANG2.0/d3yoon/1Dtest/3rdTestInvs/ArtNet_place.def
} else {
  puts "Invalid testcase"
  exit
}

# report_pack_hpwl
detailed_placement
# report_pack_hpwl
# set movecount $::env(MOVE)
# set sa_params_file "saParam_${movecount}M.json"
setSAParams -json_file "./setSAParam.json" 

opt_sa_1d
check_placement

detailed_placement

write_def ArtNet_dpo.def

set stamp2 [clock seconds]

puts "\[INFO\] Running time:   [expr $stamp2 - $stamp1] seconds"

exit

