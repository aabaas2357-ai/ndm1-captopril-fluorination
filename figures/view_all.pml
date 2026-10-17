reinitialize
load target_5ZJ2.pdb, prot
hide everything, prot
show cartoon, prot
color gray80, prot
set cartoon_transparency, 0.5, prot
select zn, prot and elem Zn
show spheres, zn
color orange, zn
set sphere_scale, 0.45, zn
select key, prot and resi 67+70+73+120+122+124+189+208+211+220+250
show sticks, key
color palecyan, key and elem C
set stick_radius, 0.14
label key and name CA, resn+resi
set label_size, -0.6
bg_color white
load h1_tetrazole_local.sdf, h1_tetrazole
show sticks, h1_tetrazole
color yellow, h1_tetrazole and elem C
set stick_radius, 0.2, h1_tetrazole
dist h1_tetrazole_SZn, (h1_tetrazole and elem S), zn, 3.5
hide labels, h1_tetrazole_SZn
load h2_isonip_tz_local.sdf, h2_isonip_tz
show sticks, h2_isonip_tz
color salmon, h2_isonip_tz and elem C
set stick_radius, 0.2, h2_isonip_tz
dist h2_isonip_tz_SZn, (h2_isonip_tz and elem S), zn, 3.5
hide labels, h2_isonip_tz_SZn
load h3_thioPro_local.sdf, h3_thioPro
show sticks, h3_thioPro
color orange, h3_thioPro and elem C
set stick_radius, 0.2, h3_thioPro
dist h3_thioPro_SZn, (h3_thioPro and elem S), zn, 3.5
hide labels, h3_thioPro_SZn
load h4_indoline_local.sdf, h4_indoline
show sticks, h4_indoline
color slate, h4_indoline and elem C
set stick_radius, 0.2, h4_indoline
dist h4_indoline_SZn, (h4_indoline and elem S), zn, 3.5
hide labels, h4_indoline_SZn
load h5_thiazole_local.sdf, h5_thiazole
show sticks, h5_thiazole
color violet, h5_thiazole and elem C
set stick_radius, 0.2, h5_thiazole
dist h5_thiazole_SZn, (h5_thiazole and elem S), zn, 3.5
hide labels, h5_thiazole_SZn
load ref_CF2_44diF_local.sdf, ref_CF2_44diF
show sticks, ref_CF2_44diF
color magenta, ref_CF2_44diF and elem C
set stick_radius, 0.2, ref_CF2_44diF
dist ref_CF2_44diF_SZn, (ref_CF2_44diF and elem S), zn, 3.5
hide labels, ref_CF2_44diF_SZn
load ref_CF3_44diF_local.sdf, ref_CF3_44diF
show sticks, ref_CF3_44diF
color teal, ref_CF3_44diF and elem C
set stick_radius, 0.2, ref_CF3_44diF
dist ref_CF3_44diF_SZn, (ref_CF3_44diF and elem S), zn, 3.5
hide labels, ref_CF3_44diF_SZn
load ref_CF3chain_local.sdf, ref_CF3chain
show sticks, ref_CF3chain
color cyan, ref_CF3chain and elem C
set stick_radius, 0.2, ref_CF3chain
dist ref_CF3chain_SZn, (ref_CF3chain and elem S), zn, 3.5
hide labels, ref_CF3chain_SZn
load ref_Dcap_local.sdf, ref_Dcap
show sticks, ref_Dcap
color green, ref_Dcap and elem C
set stick_radius, 0.2, ref_Dcap
dist ref_Dcap_SZn, (ref_Dcap and elem S), zn, 3.5
hide labels, ref_Dcap_SZn
load val_Lcap_local.sdf, val_Lcap
show sticks, val_Lcap
color lime, val_Lcap and elem C
set stick_radius, 0.2, val_Lcap
dist val_Lcap_SZn, (val_Lcap and elem S), zn, 3.5
hide labels, val_Lcap_SZn
load val_RR_local.sdf, val_RR
show sticks, val_RR
color hotpink, val_RR and elem C
set stick_radius, 0.2, val_RR
dist val_RR_SZn, (val_RR and elem S), zn, 3.5
hide labels, val_RR_SZn
load val_RS_local.sdf, val_RS
show sticks, val_RS
color wheat, val_RS and elem C
set stick_radius, 0.2, val_RS
dist val_RS_SZn, (val_RS and elem S), zn, 3.5
hide labels, val_RS_SZn
load val_SR_local.sdf, val_SR
show sticks, val_SR
color yellow, val_SR and elem C
set stick_radius, 0.2, val_SR
dist val_SR_SZn, (val_SR and elem S), zn, 3.5
hide labels, val_SR_SZn
load val_SS_local.sdf, val_SS
show sticks, val_SS
color salmon, val_SS and elem C
set stick_radius, 0.2, val_SS
dist val_SS_SZn, (val_SS and elem S), zn, 3.5
hide labels, val_SS_SZn
disable all
enable prot
enable zn
enable key
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable h1_tetrazole
enable h1_tetrazole_SZn
zoom (h1_tetrazole or zn), 9
scene h1_tetrazole, store, "h1_tetrazole | local -3.938 | S-Zn 4.45 | contact 3.02"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable h2_isonip_tz
enable h2_isonip_tz_SZn
zoom (h2_isonip_tz or zn), 9
scene h2_isonip_tz, store, "h2_isonip_tz | local -5.259 | S-Zn 5.79 | contact 3.23"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable h3_thioPro
enable h3_thioPro_SZn
zoom (h3_thioPro or zn), 9
scene h3_thioPro, store, "h3_thioPro | local -3.361 | S-Zn 2.43 | contact 2.72"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable h4_indoline
enable h4_indoline_SZn
zoom (h4_indoline or zn), 9
scene h4_indoline, store, "h4_indoline | local -2.127 | S-Zn 14.19 | contact 3.05"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable h5_thiazole
enable h5_thiazole_SZn
zoom (h5_thiazole or zn), 9
scene h5_thiazole, store, "h5_thiazole | local -2.319 | S-Zn 5.32 | contact 2.9"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable ref_CF2_44diF
enable ref_CF2_44diF_SZn
zoom (ref_CF2_44diF or zn), 9
scene ref_CF2_44diF, store, "ref_CF2_44diF | local -3.268 | S-Zn 2.25 | contact 2.97"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable ref_CF3_44diF
enable ref_CF3_44diF_SZn
zoom (ref_CF3_44diF or zn), 9
scene ref_CF3_44diF, store, "ref_CF3_44diF | local -4.572 | S-Zn 2.89 | contact 2.92"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable ref_CF3chain
enable ref_CF3chain_SZn
zoom (ref_CF3chain or zn), 9
scene ref_CF3chain, store, "ref_CF3chain | local -3.737 | S-Zn 2.83 | contact 2.93"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable ref_Dcap
enable ref_Dcap_SZn
zoom (ref_Dcap or zn), 9
scene ref_Dcap, store, "ref_Dcap | local -3.879 | S-Zn 2.63 | contact 2.82"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable val_Lcap
enable val_Lcap_SZn
zoom (val_Lcap or zn), 9
scene val_Lcap, store, "val_Lcap | local -3.985 | S-Zn 2.51 | contact 2.95"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable val_RR
enable val_RR_SZn
zoom (val_RR or zn), 9
scene val_RR, store, "val_RR | local -4.273 | S-Zn 2.9 | contact 3.07"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable val_RS
enable val_RS_SZn
zoom (val_RS or zn), 9
scene val_RS, store, "val_RS | local -4.322 | S-Zn 3.38 | contact 2.98"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable val_SR
enable val_SR_SZn
zoom (val_SR or zn), 9
scene val_SR, store, "val_SR | local -3.353 | S-Zn 2.85 | contact 2.82"
disable h1_tetrazole h1_tetrazole_SZn h2_isonip_tz h2_isonip_tz_SZn h3_thioPro h3_thioPro_SZn h4_indoline h4_indoline_SZn h5_thiazole h5_thiazole_SZn ref_CF2_44diF ref_CF2_44diF_SZn ref_CF3_44diF ref_CF3_44diF_SZn ref_CF3chain ref_CF3chain_SZn ref_Dcap ref_Dcap_SZn val_Lcap val_Lcap_SZn val_RR val_RR_SZn val_RS val_RS_SZn val_SR val_SR_SZn val_SS val_SS_SZn
enable ref_Dcap
enable ref_Dcap_SZn
enable val_SS
enable val_SS_SZn
zoom (val_SS or zn), 9
scene val_SS, store, "val_SS | local -3.569 | S-Zn 2.4 | contact 3.06"
set scene_buttons, 1
scene ref_Dcap, recall
