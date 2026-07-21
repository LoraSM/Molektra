# Copyright (C) 2026 MOLEKTRA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
import sys
from PyQt5.QtGui import QSurfaceFormat
from PyQt5.QtWidgets import QApplication
from gui.main_window import MainWindow
from PyQt5.QtOpenGL import QGLFormat

def main():
    fmt = QSurfaceFormat()
    fmt.setVersion(2,1)
    fmt.setProfile(QSurfaceFormat.NoProfile)
    fmt.setSwapBehavior(QSurfaceFormat.DoubleBuffer)
    fmt.setDepthBufferSize(24)
    fmt.setStencilBufferSize(8)
    QSurfaceFormat.setDefaultFormat(fmt)
    glfmt = QGLFormat()
    glfmt.setDoubleBuffer(True)
    glfmt.setDepthBufferSize(24)
    glfmt.setStencil(True)
    glfmt.setStencilBufferSize(8)
    glfmt.setVersion(2,1)
    glfmt.setProfile(QGLFormat.NoProfile)
    QGLFormat.setDefaultFormat(glfmt)
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
