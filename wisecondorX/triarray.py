#!/usr/bin/env python

import numpy as np


class TriArr:

    def __init__(self, new_edge):
        self.edge = new_edge
        size = int((self.edge * self.edge) / 2. + self.edge / 2.)
        self.data_array = np.zeros(size)
        self.add_index = 0

    def addValue(self, value):
        self.data_array[self.add_index] = value
        self.add_index += 1

    def setValue(self, x, y, value):
        self.data_array[x * self.edge - sum(range(x)) + y - x] = value

    def getValue(self, x, y):
        return self.data_array[x * self.edge - sum(range(x)) + y - x]

    def getSubTriangle(self, start, end):
        sub = TriArr(end - start)
        subAdd = sub.addValue
        gval = self.getValue
        for x in range(start, end):
            for y in range(x, end):
                subAdd(gval(x, y))
        return sub

    def toString(self):
        outString = 'derp'
        for x in range(self.edge):
            pass
        return outString

    def linTo2D(self, y):
        curEdge = self.edge
        while y >= curEdge:
            y -= curEdge
            curEdge -= 1
        return self.edge - curEdge, y + self.edge - curEdge
        y = 0
        while z >= int((y * y) / 2. + y / 2.):
            y += 1
        y -= 1
        x = z - int((y * y) / 2. + y / 2.)
        return x, y

    def segmentTri(self, threshold, min_search=3):
        myResult = []

        champPos = self.data_array.argmax()
        champVal = self.data_array[champPos]

        botPos = self.data_array.argmin()
        botVal = self.data_array[botPos]

        if abs(botVal) > champVal:
            champVal = botVal
            champPos = botPos

        if abs(champVal) < threshold:
            return myResult

        x, y = self.linTo2D(champPos)
        if x > min_search:
            myResult.extend(self.getSubTriangle(0, x).segmentTri(threshold, min_search))
        myResult.append((champVal, (x, y)))
        if y + 1 < self.edge - min_search:
            rightEnd = self.getSubTriangle(y + 1, self.edge).segmentTri(threshold, min_search)
            rightEnd = [(z[0], (z[1][0] + y + 1, z[1][1] + y + 1)) for z in rightEnd]
            myResult.extend(rightEnd)

        return myResult
