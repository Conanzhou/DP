class menuFileAction():
    def setSystemNameTrigger(self):
        global drawMode
        myShot = self.shot.value()
        if self.mode == 0 or self.mode==1:
            DPI.setSystemName(myShot)
        elif self.mode ==2:
            pass
        elif self.mode==3:
            pass

    def saveDataTrigger(self):   # save data for later analysis
        myShot = self.shot.value()
        path = Path.cwd() / "data" / (str(myShot) + '.mat')
        dataFile,_ =QFileDialog.getSaveFileName(self,'OpenFile',path,"matlab files (*.mat)")
        myCurves=self.curveData

        sio.savemat(dataFile, {'ChnlName': myCurves.ChnlName, 'X': myCurves.X, 'Y': myCurves.Y})
        # sio.savemat(dataFile, {'ChnlName': myCurves.ChnlName, 'X': myCurves.X, 'Y': myCurves.Y, 'myCurves': myCurves})
    
    def loadDataTrigger(self):   # load data for analysis

        path = Path.cwd() / "data"
        dataFile,_ =QFileDialog.getOpenFileName(self, 'OpenFile', path, "matlab files (*.mat)")
        f=sio.loadmat(dataFile)
        ChnlName=f['ChnlName']   #
        X=f['X']   #
        Y=f['Y']   #
        # myCurves=f['myCurves']

        self.Curves.clear()
        self.curveData = CurveData('Songxm')  # initializing the class

        for n in range(len(ChnlName)):
            channelName=ChnlName[n]
            self.Curves.addItem(channelName)
            x=X[n]
            y=Y[n]
            self.setCurveData(x, y, n, channelName, channelName)

    def saveConfigTrigger(self):  # save channel configuration for later use

        path = Path.cwd() / "configuration"
        configFile,_ =QFileDialog.getSaveFileName(self,'OpenFile',path,"matlab files (*.mat)")
        myCurves=self.curveData
        myCurves.x=[]
        myCurves.y=[]
        sio.savemat(configFile, {'ChnlName': myCurves.ChnlName})
    

    def loadConfigTrigger(self):
        path = Path.cwd() / "configuration"
        configFile,_ =QFileDialog.getOpenFileName(self,'OpenFile',path,"matlab files (*.mat)")
        f=sio.loadmat(configFile)
        Chnl=f['ChnlName']   # np.array

        self.hD = selectChnl(Chnl)
        self.hD.setWindowModality(Qt.ApplicationModal)
        self.hD.show()
        self.hD.signalChnl.connect(self.receiveChnl)

    
    def drawModeTrigger(self):
        global drawMode
        number, ok = QInputDialog.getInt(self, "action", "0=mat,1=pyqtgrapg", drawMode, 0, 1, 1)
        if ok:
            drawMode = number

        self.draw.clicked.disconnect()

        if drawMode == 0:  # matlab mode
            self.draw.clicked.connect(self.drawClickedMat)
        elif drawMode == 1:  # pyqtgraph mode
            self.draw.clicked.connect(self.drawClicked)

class menumachine:
    def setMachineTrigger():
        global setMachineMode
        setMachineMode=1
        self.Warning.clear()
        self.Warning.append('please select the data source')
        machineNames=['hl2a', 'localdas', 'east', 'exl50', 'hl2m']
        self.Files.clear()
        # show the machine in Files for select
        for m in machineNames:
            self.Files.addItem(m)

    def ChnlConfigTrigger():  # modify the configuration
        self.hD = modifyConfig(self.curveData)
        self.hD.setWindowModality(Qt.ApplicationModal)
        self.hD.show()
        self.hD.signalCurve.connect(self.receiveCurveData)
    
    def scanTreeTrigger():  # scan tree for channels
        myShot = self.shot.value()

        if self.mode == 0 or self.mode == 1:
            pass
        elif self.mode == 2:
            DPI.sct('exl50', myShot)
        elif self.mode == 3:
            DPI.sct('east', myShot)

    def scanTreeADTrigger():  # scan tree for channels
        myShot = self.shot.value()

        if self.mode == 0 or self.mode == 1:
            pass
        elif self.mode == 2:
            DPI.sctAD('exl50', myShot)
        elif self.mode == 3:
            DPI.sctAD('east', myShot)

    def defaultChnlTrigger():  # modify the configuration
        self.defaultChnlSaveLoad()  # save

class menuHelp:
    def BasicsTrigger():  # how to use
        self.Warning.append('---------------------------------------------------------------------------------' +
                            '\n[BASICS]\n' +
                            'Choosing a machine: click Machine and select setMachine from the drop down bar. Choose '
                            + 'desired machine listed in the Files module\n\n' +
                            'Searching for channel patterns: input regular expressions into search bar labeled '
                            + 'Search Chnl. Results will appear in the Browser module\n\n' +
                            'Adding channels to Curves: click Browser to change it to Browser+Add. You can now '
                            + 'select multiple channels which will be listed in the Curves module, and can be drawn, '
                            + 'updated, exported, etc\n\n' +
                            'Exporting/Importing channels: click File and select saveConfig/loadConfig\n\n' +
                            'For more information, click Help\n' +
                            '---------------------------------------------------------------------------------')
    def ButtonsTrigger():
        self.Warning.append('---------------------------------------------------------------------------------' +
                            '\n[BUTTONS]\n' +
                            'shotOK: prepares work for data mining for shot\n' +
                            'update: reloads data based on current settings\n' +
                            'draw: displays multiple data side by side when in Browser+Add mode\n' +
                            'clearAll: clears entire interface\n' +
                            '---------------------------------------------------------------------------------')
    def DropDownMenuTrigger():
        self.Warning.append('---------------------------------------------------------------------------------' +
                            '\n[DROP DOWN MENU]\n' +
                            'setSystemName: [INFO]\n' +
                            'saveConfig: saves current settings into external file\n' +
                            'saveData: [INFO]\n' +
                            'loadConfig: loads selected external file to update current settings\n' +
                            'setMachine: allows user to select machine, listed under Files\n' +
                            'ChnlConfig: allows user to modify the channel information\n' +
                            'defaultChnl: [INFO]\n' +
                            '---------------------------------------------------------------------------------')
