mapkey gencam @MAPKEY_NAMEGenerate 3D shape from cam profile data in \
mapkey(continued) camProfile.pts in working directory;@MAPKEY_LABELGenerate Cam;\
mapkey(continued) ~ Command `ProCmdModelNew` ;\
mapkey(continued) ~ Activate `new` `chk_use_default_template` 0;~ Activate `new` `OK`;\
mapkey(continued) ~ FocusOut `new_file_opts` `inp_template_name`;\
mapkey(continued) ~ Activate `new_file_opts` `psh_browse`;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `DLG_PREVIEW_POST` `file_open`;\
mapkey(continued) ~ Activate `file_open` `pb_favorites__FAV_19_`;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `PREVIEW_POPUP_TIMER` \
mapkey(continued) `file_open:Ph_list.Filelist:<NULL>`;\
mapkey(continued) ~ Select `file_open` `Ph_list.Filelist` 1 `camtemplate.prt`;\
mapkey(continued) ~ Command `ProFileSelPushOpen_Standard@context_dlg_open_cmd` ;\
mapkey(continued) ~ FocusIn `new_file_opts` `psh_browse`;~ Activate `new_file_opts` `psh_ok`;\
mapkey(continued) ~ Command `ProCmdDatumPointOffsetCSys` ;\
mapkey(continued) ~ Select `main_dlg_cur` `PHTLeft.AssyTree` 1 `node7`;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `PREVIEW_POPUP_TIMER` \
mapkey(continued) `main_dlg_w4:PHTLeft.AssyTree:<NULL>`;\
mapkey(continued) ~ Activate `Odui_Dlg_00` `t1.ImportBtn`;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `DLG_PREVIEW_POST` `file_open`;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `PREVIEW_POPUP_TIMER` \
mapkey(continued) `file_open:Ph_list.Filelist:<NULL>`;\
mapkey(continued) ~ Select `file_open` `Ph_list.Filelist` 1 `camDirection.pts`;\
mapkey(continued) ~ Command `ProFileSelPushOpen_Import@context_dlg_open_cmd` ;\
mapkey(continued) ~ Activate `Odui_Dlg_00` `stdbtn_1`;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `PREVIEW_POPUP_TIMER` \
mapkey(continued) `main_dlg_w4:PHTLeft.AssyTree:<NULL>`;\
mapkey(continued) ~ Command `ProCmdDatumPointOffsetCSys` ;\
mapkey(continued) ~ Select `main_dlg_cur` `PHTLeft.AssyTree` 1 `node7`;\
mapkey(continued) ~ Activate `Odui_Dlg_00` `t1.ImportBtn`;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `DLG_PREVIEW_POST` `file_open`;\
mapkey(continued) ~ Select `file_open` `Ph_list.Filelist` 1 `camProfile.pts`;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `PREVIEW_POPUP_TIMER` \
mapkey(continued) `file_open:Ph_list.Filelist:<NULL>`;\
mapkey(continued) ~ Command `ProFileSelPushOpen_Import@context_dlg_open_cmd` ;\
mapkey(continued) ~ Activate `Odui_Dlg_00` `stdbtn_1`;~ Command `ProCmdFtCrvTp` ;\
mapkey(continued) ~ Trail `UI Desktop` `UI Desktop` `PREVIEW_POPUP_TIMER` \
mapkey(continued) `main_dlg_w4:PHTLeft.AssyTree:<NULL>`;\
mapkey(continued) ~ Activate `main_dlg_cur` `dashInst0.Done`;~ Command `ProCmdFtExtrude` ;\
mapkey(continued) ~ Command `ProCmdEnvDtmDisp`  1;@PAUSE_FOR_SCREEN_PICK;\
mapkey(continued) ~ Command `ProCmdSketOffset`  1;@PAUSE_FOR_SCREEN_PICK;0;\
mapkey(continued) ~ Activate `useedge` `CloseButton`;~ Command `ProCmdSketDone`;