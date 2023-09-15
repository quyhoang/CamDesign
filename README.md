# Cam Design Assistant


## Introduction

A program supporting mechanical designers in the seamless development of planar oscillating and translating cam systems.

## Table of Contents

1. [Introduction](#introduction)
2. [Table of Contents](#table-of-contents)
3. [Features](#features)
4. [Functionality](#functionality)
5. [Usage](#usage)
6. [Contributing](#contributing)
7. [License](#license)

## Features

- Simple syntax
- Easy to use (the program is written with Matlab OOP)
- Free (can be used with Matlab online)

## Functionality

CamDesign helps you to investigate and design the following:

- Cam follower and load:
  - Displacement
  - Velocity
  - Acceleration

- Radius of curvature 
- Pressure angle
- Pitch curve
- Cam profile
- Cam machining process
- Cam motion
- Spring force
- Motor torque

## Usage

### How to initialize a cam object

First, create 2 configuration files:

**&lt;configname&gt;_transition.txt**

This file contains information about the trajectory of load or cam follower.

**&lt;configname&gt;_parameter.txt**

This file contains all the parameters needed to define a cam system.

Replace &lt;configname&gt; with the actual name of your configuration.
You can find examples of configuration files in source code folder.

After you have 2 configuration files with same configuration name, to create a cam object, just initiate the object using either **ocam** or **tcam**.

ocam is for oscillating cam. tcam is for translating cam.

Example:

&gt;&gt;mycam = ocam('mycamname');

The command above will looking for mycamname_transition.txt and mycamname_parameter.txt and create a cam object using the information in those two files.



### How to investigate the cam object



#### camdata class object

The following attributes of the cam object are instance of camdata class.
They have built-in methods to show, plot, and export their content to xlsx files:

- displacement
- load_displacement
- velocity
- acceleration
- curvature
- pressureAngle
- motorTorque



To view their name, unit, and data, just type their name to the Matlab command console: 
**obj**

To view their data, use the attribute data: 
**obj.data**

To view the plot of their data vs rotation angle of the cam, use show method: 
**obj.show**

To export the data to xlsx file, use the excel method: 
**obj.export**

To summarize, for example, if you work with the cam object named cam and the attribute velocity.

| Command Console Syntax  | Function  | 
| :------------: |:---------------:| 
| cam.velocity | show name, unit, and data |  
| cam.velocity.data      | show numerical data        |        
| cam.velocity.show     | show plot   |    
| cam.velocity.export     | export to xlsx   |    


Because some plot are frequenly viewed, there are shortcuts to quickly create their plot.

For example, if the cam object name is cam.

| Data to plot  | Standard method | Shortcut |
| :------------: |:---------------:| :-------------:|
| load_displacement | cam.load_displacement.show |   cam.s  |
| velocity      | cam.velocity.show        |          cam.v  |
| acceleration     | cam.acceleration.show   |          cam.a    |
| curvature     | cam.curvature.show   |          cam.c   |
| pressureAngle     | cam.pressureAngle.show   |          cam.pressure    |
| motorTorque     | cam.motorTorque.show   |          cam.torque    |

Example:

**&gt;&gt;cam.s**

![displacement](https://github.com/quyhoang/CamDesign/assets/14304980/070d9390-c8ca-40ec-881e-8279120a2252)



#### other attributes

Some of frequenly used attributes

| Attribute name  | Explanation  | 
| :------------: |:---------------:| 
| minK | minimum acceptable spring constant |  
| fSpring      | spring force        |        
| fCut     | additional force needed for cutting or pressing   |      




#### Methods

**Class methods**

| Methods  | Explanation  | 
| :------------: |:---------------:| 
| camCheck | Quick check whether the cam is valid |  
| profile      | Show cam profile        |        
| export     | Export profile data as points. These data can be used to directly generate 3D cam model using CreoAutomation (another software written by me)  |    
| machining     | Show machining process   | 
| spring     | Plot spring force   | 
| animation     | Animating cam and follower motion   | 

Example:

**camCheck**

&gt;&gt; cam.camCheck
 
最大速度: 188.2804 mm/s

最大加速: 15.7934 m/s2

最大圧力角: -9.422 °

最大モータートルク: 2.056 N/m

最小曲率半径: 35.2866 mm

最小許容バネ定数: 0.26751 N/mm

**machining**

&gt;&gt; cam.machining

![machining](https://github.com/quyhoang/CamDesign/assets/14304980/37c860e6-1c02-4926-b0b2-8c1b002593fd)

**animation**

&gt;&gt; cam.animation

![animation](https://github.com/quyhoang/CamDesign/assets/14304980/2901a0c3-8ecc-4e01-9e58-96f2e816ddf8)


**Independent methods**

| Methods  | Explanation  | 
| :------------: |:---------------:| 
| combine | Calculate sum and plot all data over the period of 360 degree |  
| excel | Combine and export all data field of input arguments to xlsx |  

Example (number of inputs is arbitrary):

&gt;&gt; combine(cam1.torque.data,cam2.torque.data,cam3.torque.data)

![Combine](https://github.com/quyhoang/CamDesign/assets/14304980/74afc908-dc0e-45c5-b250-a068f1aafb32)

&gt;&gt; excel('outputfile',cam1.motorTorque,cam1.velocity,cam1.displacement)

<p align="center">
  <img src="https://github.com/quyhoang/CamDesign/assets/14304980/6075e723-1075-4ba2-8d57-ce95c1ded13b" alt="Excel" width="400" />
</p>



## Contributing

To contribute to this project, please contact me first.

## License

This project is licensed under the terms of the GNU General Public License v3.0.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)



Thank you for visiting this project. Feel free to contribute and make this project better.


---

# カムデザインアシスタントプログラム

(機械翻訳)

## はじめに

平面振動カムおよび移動カムの設計を支援するプログラム。

## 目次

1. [はじめに](#はじめに)
2. [目次](#目次)
3. [機能](#機能)
4. [機能性](#機能性)
5. [使い方](#使い方)
6. [貢献](#貢献)
7. [ライセンス](#ライセンス)

## 機能

- シンプルな構文
- 使いやすい（プログラムはMatlab OOPで書かれています）
- 無料（Matlabオンラインで使用可能）

## 機能性

CamDesignは以下の調査と設計を支援します：

- カムフォロワーと負荷：
  - 変位
  - 速度
  - 加速度

- 曲率半径
- 圧力角
- ピッチ曲線
- カムプロフィール
- カム加工プロセス
- カムモーション
- バネ力
- モータートルク

## 使い方

### カムオブジェクトの初期化方法

まず、2つの構成ファイルを作成します：

**&lt;configname&gt;_transition.txt**

このファイルには、負荷またはカムフォロワーの軌跡に関する情報が含まれています。

**&lt;configname&gt;_parameter.txt**

このファイルには、カムを作成するために必要なすべてのパラメータが含まれています。

&lt;configname&gt;を実際の構成名に置き換えます。
ソースコードフォルダに構成ファイルの例があります。

同じ構成名で2つの構成ファイルが作成された後、カムオブジェクトを作成するには、**ocam**または**tcam**を使用してオブジェクトを開始するだけです。

ocamは振動カム用です。tcamは移動カム用です。

例：

&gt;&gt;mycam = ocam('mycamname');

上記のコマンドは、mycamname_transition.txtとmycamname_parameter.txtを検索し、それら2つのファイルの情報を使用してカムオブジェクトを作成します。

.

### カムオブジェクトの調査方法

#### camdataクラスオブジェクト

カムオブジェクトの次の属性はcamdataクラスのインスタンスです。
xlsxファイルに内容を表示、プロット、エクスポートする組み込みメソッドがあります：

- 変位
- 負荷変位
- 速度
- 加速度
- 曲率
- 圧力角
- モータートルク

Matlabコマンドコンソールに名前、単位、データを表示するには、その名前を入力するだけです：

**obj**

データを表示するには、属性データを使用します：

**obj.data**

カムの回転角に対するデータのプロットを表示するには、showメソッドを使用します：

**obj.show**

データをxlsxファイルにエクスポートするには、excelメソッドを使用します：
**obj.export**

例として、カムオブジェクト名がcamで、属性が速度の場合をまとめます。

| コマンドコンソールの構文  | 機能  | 
| :------------: |:---------------:| 
| cam.velocity | 名前、単位、データを表示 |  
| cam.velocity.data      | 数値データを表示        |        
| cam.velocity.show     | プロットを表示   |    
| cam.velocity.export     | xlsxにエクスポート   |    

いくつかのプロットは頻繁に表示されるため、それらのプロットを素早く作成するショートカットがあります。

例えば、カムオブジェクト名がcamの場合。

| データをプロット  | 標準方法 | ショートカット |
| :------------: |:---------------:| :-------------:|
| load_displacement | cam.load_displacement.show |   cam.s  |
| velocity      | cam.velocity.show        |          cam.v  |
| acceleration     | cam.acceleration.show   |          cam.a    |
| curvature     | cam.curvature.show   |          cam.c   |
| pressureAngle     | cam.pressureAngle.show   |          cam.pressure    |
| motorTorque     | cam.motorTorque.show   |          cam.torque    |

例:

**&gt;&gt;cam.s**

![displacement](https://github.com/quyhoang/CamDesign/assets/14304980/070d9390-c8ca-40ec-881e-8279120a2252)


#### 他の属性

頻繁に使用されるいくつかの属性

| 属性名  | 説明  | 
| :------------: |:---------------:| 
| minK | 最小許容バネ定数 |  
| fSpring      | バネ力        |        
| fCut     | 切削やプレスに必要な追加力   |      



#### メソッド

**クラスメソッド**

| メソッド  | 説明  | 
| :------------: |:---------------:| 
| camCheck | カムが有効かどうかを素早く確認 |  
| profile      | カムプロファイルを表示        |        
| export     | プロファイルデータを点としてエクスポート。これらのデータは、CreoAutomation（私が書いた別のソフトウェア）を使用して直接3Dカムモデルを生成するために使用できます  |    
| machining     | 加工プロセスを表示   | 
| spring     | バネ力をプロット   | 
| animation     | カムとフォロワーの動きをアニメーション化   | 

例：

**camCheck**

&gt;&gt; cam.camCheck
 
最大速度: 188.2804 mm/s

最大加速: 15.7934 m/s2

最大圧力角: -9.422 °

最大モータートルク: 2.056 N/m

最小曲率半径: 35.2866 mm

最小許容バネ定数: 0.26751 N/mm

**machining**

&gt;&gt; cam.machining

![machining](https://github.com/quyhoang/CamDesign/assets/14304980/37c860e6-1c02-4926-b0b2-8c1b002593fd)

**animation**

&gt;&gt; cam.animation

![animation](https://github.com/quyhoang/CamDesign/assets/14304980/2901a0c3-8ecc-4e01-9e58-96f2e816ddf8)


**独立した方法**

| メソッド  | 説明  | 
| :------------: |:---------------:| 
| combine | 360度の期間にわたるすべてのデータの合計とプロットを計算 |  
| excel | 入力引数のすべてのデータフィールドを組み合わせてxlsxにエクスポート |  

例（入力の数は任意）：

&gt;&gt; combine(cam1.torque.data,cam2.torque.data,cam3.torque.data)

![Combine](https://github.com/quyhoang/CamDesign/assets/14304980/74afc908-dc0e-45c5-b250-a068f1aafb32)

&gt;&gt; excel('outputfile',cam1.motorTorque,cam1.velocity,cam1.displacement)

<p align="center">
  <img src="https://github.com/quyhoang/CamDesign/assets/14304980/6075e723-1075-4ba2-8d57-ce95c1ded13b" alt="Excel" width="400" />
</p>



## 貢献

このプロジェクトに貢献するには、最初に私に連絡してください。

## ライセンス

このプロジェクトはGNU General Public License v3.0の条件の下でライセンスされています。

[![ライセンス: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

このプロジェクトを訪れていただきありがとうございます。気軽に貢献して、このプロジェクトをより良くしてください。


---
