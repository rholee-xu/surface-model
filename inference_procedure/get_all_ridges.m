% get all outline

% xy ridge caulonema
bxy_220830_c3 = get_ridge('Raw Ridge Files/220830_xyridge_before_c3_s72.csv',72);
axy_220830_c3 = get_ridge('Raw Ridge Files/220830_xyridge_after_c3_s72.csv',72);

bxy_220825_c4 = get_ridge('Raw Ridge Files/220825_xyridge_before_c4.csv',73);
axy_220825_c4 = get_ridge('Raw Ridge Files/220825_xyridge_after_c4.csv',58);

bxy_221006_c1 = get_ridge('Raw Ridge Files/221006_xyridge_before_c1_s56.csv',56);
axy_221006_c1 = get_ridge('Raw Ridge Files/221006_xyridge_after_c1_s57.csv',57);

bxy_220825_c3 = get_ridge('Raw Ridge Files/220825_xyridge_before_c3_s70.csv',70);
axy_220825_c3 = get_ridge('Raw Ridge Files/220825_xyridge_after_c3_s70.csv',70);

bxy_221102_c3 = get_ridge('Raw Ridge Files/221102_xyridge_before_c3_s71.csv',71);
axy_221102_c3 = get_ridge('Raw Ridge Files/221102_xyridge_after_c3_s68.csv',68);

bxy_221102_c4 = get_ridge('Raw Ridge Files/221102_xyridge_before_c4_s76.csv',76);
axy_221102_c4 = get_ridge('Raw Ridge Files/221102_xyridge_after_c4_s74.csv',74);

compareplot(bxy_220830_c3,bxy_220825_c3); hold on; compareplot(bxy_221006_c1,bxy_220825_c4);
hold on; compareplot(bxy_221102_c4,bxy_221102_c3);

% xy ridge Chloronema
bxy_221003_c3 = get_ridge('Raw Ridge Files/221003_xyridge_before_c3_s60.csv',60);
axy_221003_c3 = get_ridge('Raw Ridge Files/221003_xyridge_after_c3_s50.csv',50);

bxy_221006_c3 = get_ridge('Raw Ridge Files/221006_xyridge_before_c3_s92.csv',92);
axy_221006_c3 = get_ridge('Raw Ridge Files/221006_xyridge_after_c3_s67.csv',67);

bxy_221007_c4 = get_ridge('Raw Ridge Files/221007_xyridge_before_c4_s64.csv',64);
axy_221007_c4 = get_ridge('Raw Ridge Files/221007_xyridge_after_c4_s43.csv',43);

bxy_220825_c2 = get_ridge('Raw Ridge Files/220825_xyridge_before_c2_s70.csv',70);
axy_220825_c2 = get_ridge('Raw Ridge Files/220825_xyridge_after_c2_s77.csv',77);

bxy_220829_c1 = get_ridge('Raw Ridge Files/220829_xyridge_before_c1_s95.csv',95);
axy_220829_c1 = get_ridge('Raw Ridge Files/220829_xyridge_after_c1_s93.csv',93);

bxy_221027_c4 = get_ridge('Raw Ridge Files/221027_xyridge_before_c4_s78.csv',78);
axy_221027_c4 = get_ridge('Raw Ridge Files/221027_xyridge_after_c4_s76.csv',76);

bxy_221027_c5 = get_ridge('Raw Ridge Files/221027_xyridge_before_c5_s58.csv',58);
axy_221027_c5 = get_ridge('Raw Ridge Files/221027_xyridge_after_c5_s87.csv',87);

compareplot(bxy_221003_c3,bxy_221006_c3); hold on; compareplot(bxy_221007_c4,bxy_220825_c2);
hold on; compareplot(bxy_220829_c1,bxy_221027_c4); hold on; quickplot(bxy_221027_c5);

% xz ridge caulonema
bxz_220830_c3 = get_ridge('Raw Ridge Files/220830_xzridge_before_c3_s81.csv',81);
axz_220830_c3 = get_ridge('Raw Ridge Files/220830_xzridge_after_c3_s86.csv',86);
%bxz_220830_c3 = bxz_220830_c3([1,3,2],:); axz_220830_c3 = axz_220830_c3([1,3,2],:);

bxz_220825_c4 = get_ridge('Raw Ridge Files/220825_xzridge_before_c4.csv',72);
axz_220825_c4 = get_ridge('Raw Ridge Files/220825_xzridge_after_c4.csv',72);
%bxz_220825_c4 = bxz_220825_c4([1,3,2],:); axz_220825_c4 = axz_220825_c4([1,3,2],:);

bxz_221006_c1 = get_ridge('Raw Ridge Files/221006_xzridge_before_c1_s79.csv',79);
axz_221006_c1 = get_ridge('Raw Ridge Files/221006_xzridge_after_c1_s70.csv',70);
%bxz_221006_c1 = bxz_221006_c1([1,3,2],:); axz_221006_c1 = axz_221006_c1([1,3,2],:);

bxz_220825_c3 = get_ridge('Raw Ridge Files/220825_xzridge_before_c3_s72.csv',72);
axz_220825_c3 = get_ridge('Raw Ridge Files/220825_xzridge_after_c3_s67.csv',67);
%bxz_220825_c3 = bxz_220825_c3([1,3,2],:); axz_220825_c3 = axz_220825_c3([1,3,2],:);

bxz_221102_c3 = get_ridge('Raw Ridge Files/221102_xzridge_before_c3_s28.csv',28);
axz_221102_c3 = get_ridge('Raw Ridge Files/221102_xzridge_after_c3_s24.csv',24);
%bxz_221102_c3 = bxz_221102_c3([1,3,2],:); axz_221102_c3 = axz_221102_c3([1,3,2],:);

bxz_221102_c4 = get_ridge('Raw Ridge Files/221102_xzridge_before_c4_s57.csv',57);
axz_221102_c4 = get_ridge('Raw Ridge Files/221102_xzridge_after_c4_s63.csv',63);
%bxz_221102_c4 = bxz_221102_c4([1,3,2],:); axz_221102_c4 = axz_221102_c4([1,3,2],:);

% xz ridge chloronema 
bxz_221003_c3 = get_ridge('Raw Ridge Files/221003_xzridge_before_c3_s71.csv',71);
axz_221003_c3 = get_ridge('Raw Ridge Files/221003_xzridge_after_c3_s71.csv',71);
%bxz_221003_c3 = bxz_221003_c3([1,3,2],:); axz_221003_c3 = axz_221003_c3([1,3,2],:);

bxz_221006_c3 = get_ridge('Raw Ridge Files/221006_xzridge_before_c3_s76.csv',76);
axz_221006_c3 = get_ridge('Raw Ridge Files/221006_xzridge_after_c3_s78.csv',78);
%bxz_221006_c3 = bxz_221006_c3([1,3,2],:); axz_221006_c3 = axz_221006_c3([1,3,2],:);

bxz_221007_c4 = get_ridge('Raw Ridge Files/221007_xzridge_before_c4_s82.csv',82);
axz_221007_c4 = get_ridge('Raw Ridge Files/221007_xzridge_after_c4_s78.csv',78);
%bxz_221007_c4 = bxz_221007_c4([1,3,2],:); axz_221007_c4 = axz_221007_c4([1,3,2],:);

bxz_220825_c2 = get_ridge('Raw Ridge Files/220825_xzridge_before_c2_s83.csv',83);
axz_220825_c2 = get_ridge('Raw Ridge Files/220825_xzridge_after_c2_s84.csv',84);
%bxz_220825_c2 = bxz_220825_c2([1,3,2],:); axz_220825_c2 = axz_220825_c2([1,3,2],:);

bxz_220829_c1 = get_ridge('Raw Ridge Files/220829_xzridge_before_c1_s75.csv',75);
axz_220829_c1 = get_ridge('Raw Ridge Files/220829_xzridge_after_c1_s96.csv',96);
%bxz_220829_c1 = bxz_220829_c1([1,3,2],:); axz_220829_c1 = axz_220829_c1([1,3,2],:);

bxz_221027_c4 = get_ridge('Raw Ridge Files/221027_xzridge_before_c4_s74.csv',74);
axz_221027_c4 = get_ridge('Raw Ridge Files/221027_xzridge_after_c4_s79.csv',79);
%bxz_221027_c4 = bxz_221027_c4([1,3,2],:); axz_221027_c4 = axz_221027_c4([1,3,2],:);

bxz_221027_c5 = get_ridge('Raw Ridge Files/221027_xzridge_before_c5_s61.csv',61);
axz_221027_c5 = get_ridge('Raw Ridge Files/221027_xzridge_after_c5_s63.csv',63);
%bxz_221027_c5 = bxz_221027_c5([1,3,2],:); axz_221027_c5 = axz_221027_c5([1,3,2],:);

% plotting
compareplot(bxy_220830_c3,bxz_220830_c3);
compareplot(bxy_220825_c4,bxz_220825_c4);
compareplot(bxy_221006_c1,bxz_221006_c1);
compareplot(bxy_220825_c3,bxz_220825_c3);
compareplot(bxy_221102_c3,bxz_221102_c3);
compareplot(bxy_221102_c4,bxz_221102_c4);

compareplot(bxy_221003_c3,bxz_221003_c3);
compareplot(bxy_221006_c3,bxz_221006_c3);
compareplot(bxy_221007_c4,bxz_221007_c4);
compareplot(bxy_220825_c2,bxz_220825_c2);
compareplot(bxy_220829_c1,bxz_220829_c1);
compareplot(bxy_221027_c4,bxz_221027_c4);
compareplot(bxy_221027_c5,bxz_221027_c5);