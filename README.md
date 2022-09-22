# folded
Folded: an ImageJ/Fiji toolkit to describe mammalian herbivore dentition from 2D images

**Authors**: [Oscar Sanisidro](https://scholar.google.es/citations?user=TLYIuyEAAAAJ&hl=es), [Ignacio Arganda-Carreras](https://www.ikerbasque.net/en/ignacio-arganda-carreras), and [Juan L. Cantalapiedra](https://scholar.google.es/citations?hl=es&user=19hBfQ8AAAAJ)

## Installation instructions
‘folded’ is a script of the open-source graphic analysis software [Fiji](https://fiji.sc/), which is, in turn, a distribution of [ImageJ](http://rsb.info.nih.gov/nih-image). ‘folded’ has been scripted in the BeanShell language.
### Fiji installation
You can find the latest version of Fiji in the following [link](https://imagej.net/software/fiji/) together with the installation instructions and system requirements for the current version. At the moment of this publication Fiji is supported for Windows XP, Vista, 7, 8, 10, 11, etc. (in both 32 and 64-bit versions), Mac OS X 10.8 “Mountain Lion” or later, and Linux on amd64 and x86 architectures. To make sure the latest version is installed, it is recommended to run the [Updater](https://imagej.net/plugins/updater) when the program is executed for the first time. The update option can be found in the [Help / Update...] menu.

### Installing dependencies
‘folded’ makes use of two plugins that are not included in the base Fiji distribution. These are: [‘OrientationJ’](http://bigwww.epfl.ch/demo/orientation/) and [‘MorphoLibJ’](https://imagej.net/plugins/morpholibj). To install them,

1. Go to the [Help > Update...] menu
2. Select ‘Manage update sites’
3. From the list within, select ‘BIG-EPFL’ and ‘IJPB-plugins’
4. Click on 'Close'
5. Click on ‘Apply changes’
6. Restart Fiji.

Additional information on how to install plugins in Fiji manually can be found in the following link: https://imagej.net/plugins/.

### 'folded' script installation
Once Fiji and the two auxiliary plugins are installed, the file ‘folded_.bsh’ can be downloaded from the [Github repository](https://github.com/iarganda/folded). Next, copy the file to the directory \<Fiji root\>/plugins/Scripts/Plugins/Analyze/ (you may need to create it) and restart Fiji. The script then appears as a new command under the menu [Plugins/Analyze] (in the last position of the menu). 

If you are not able to find the root directory of your Fiji installation or you don't permissions to copy the file there, do not worry. You can open the script Editor ([File > New > Script]) and load the file into the script editor by the [File > Open] menu of the editor. Then click on the “Run” button and the script will be executed.
