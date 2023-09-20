#NoEnv  ; Recommended for performance and compatibility with future AutoHotkey releases.
SendMode Input  ; Recommended for new scripts due to its superior speed and reliability.
SetWorkingDir %A_ScriptDir%  ; Ensures a consistent starting directory.
#SingleInstance force  ; force new attempt to run the program if it is already running
SetTitleMatchMode, 2
FileEncoding, UTF-8

; CAM DESIGN

::gcam:: ;generate 3D Cam Model from txt file in working directory
generatecamprofile:
generatecamprofile()
return

generatecamprofile()
{
	; set working directory to Creo working directory
	wdir := getWorkingDirectory()
	SetWorkingDir %wdir%
	FileMove, camProfile.txt, camProfile.pts, 1
	FileMove, camDirection.txt, camDirection.pts, 1

	If !FileExist("camDirection.pts")
	{
		Msgbox,48,Warning, camDirection.pts は %wdir% にありません。,20
	}
	
	; Rename file, overwrite if destination file has already existed
	If !FileExist("camProfile.pts")
	{
		Msgbox,48,エラー, camProfile.pts は %wdir% にありません。,20
		return
	}
	SetWorkingDir %A_ScriptDir%
	SendInput gencam ; mapkey in Creo to read .pts file and generate 3D shape from curve
	
	return
}