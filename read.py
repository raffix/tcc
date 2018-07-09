import ogr
import shapely
from shapely.geometry import *
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import math
import json
import sys
import time
from numpy import ones,vstack
from numpy.linalg import lstsq
from collections import namedtuple

DEG = 57.2958

class Ponto:
    def __init__(self, nodo):
        self.x = nodo[0]
        self.y = nodo[1]

# GMM
class GMM:
    def __init__(self):
        self.mapa = ''
        self.a = 5
        self.b = 2
        self.Ap = 10
        #raio direto para coordenadas (aproximadamente 50m)
        self.bearing = 0.001
        self.Ah = self.a*self.Ap
        self.Arp = self.b*self.Ap

    def square(self, x, y):
        coords={}
        coords['minLongitude'] = float(y) - float(self.bearing)
        coords['maxLongitude'] = float(y) + float(self.bearing)
        coords['minLatitude'] = float(x) - float(self.bearing)
        coords['maxLatitude'] = float(x) + float(self.bearing)
        return coords

    def getAngle(self, p0, p1, p2, p3):
        v0 = np.array(p1) - np.array(p0)
        v1 = np.array(p2) - np.array(p3)
        angle = np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
        return np.degrees(angle)

    def euclid(self, A, B):
        pa = np.array(A)
        pb = np.array(B)
        return np.linalg.norm(pb-pa)
        
    def getPointsNear(self, ponto, candidatos):
        distancias = []
        for x in candidatos:
            dist = self.euclid(ponto, x)
            distancias.append([x, dist])
        distancias = sorted(distancias, key=lambda item: item[1])
        return distancias

    def getNearFromLink(self, ponto, candidato):
        candidato["geometrySorted"] = mm.getPointsNear(ponto, candidato['geometry'])
        return candidato

    def deltaB(self, headAngle, streetAngle):
        delta = headAngle - streetAngle
        if(-180 <= delta <= 180):
            return delta
        elif(delta> 180):
            return 360-delta
        else:
            return 360+delta

    def WSh(self, cos):
        return self.Ah * cos

    def D(self, p1, p2, p3):
        d = p3.x*(p1.y - p2.y) - p3.y*(p1.x - p2.x) + (p1.x*p2.y - p2.x*p1.y)
        temp = (((p1.x-p2.x)**2 + (p1.y - p2.y)**2)**0.5)
        if(temp == 0):
            d = d/temp
        return d

    def WSpd(self, D):
        return self.Ap/D

    def WSpi(self, theta):
        return self.Ap * np.cos(theta)

    def WSrp(self, alfa):
        return self.Arp*np.cos(alfa)

    def TWS(self, WSh, WSpd, WSpi, WSrp):
        return WSh+(WSpd+WSpi)+WSrp

    def getPoint(self, coords, index):
        position = [coords.latitude, coords.longitude, coords.heading]
        if(index == 0):
            return position
        return position

    def deltaE(self, vel, theta):
    	return (1*vel)*np.sin(theta)

    def deltaN(self, vel, theta):
    	return (1*vel)*np.cos(theta)

    def E(self, vel, theta, Pe):
    	return Pe +self.deltaE(vel, theta)

    def N(self, vel, theta, Pe):
    	return Pe +self.deltaN(vel, theta)
# Mapa

def getTraces():
    with open(sys.argv[0]) as f:
        data = json.load(f)
    return data

def getMap():
    driver = ogr.GetDriverByName('OSM')
    #data = driver.Open('mapa_chapeco.osm')
    data = driver.Open('mapaChapeco.osm')
    #data = driver.Open('mapaChapeco2.osm')

    layer = data.GetLayer('points')
    features=[x for x in layer]
    lines = data.GetLayer('lines')
    lines = [x for x in lines]
    mapa = {}


    data_list=[]
    for feature in features:
        data=feature.ExportToJson(as_object=True)
        coords=data['geometry']['coordinates']
        shapely_geo=Point(coords[0],coords[1])
        name=data['properties']['name']
        highway=data['properties']['highway']
        other_tags=data['properties']['other_tags']
        if other_tags and 'amenity' in other_tags:
            feat=[x for x in other_tags.split(',') if 'amenity' in x][0]
            amenity=feat[feat.rfind('>')+2:feat.rfind('"')]
        else:
            amenity=None
        data_list.append([name,highway,amenity,shapely_geo])

    mapa['pontos'] = gpd.GeoDataFrame(data_list,columns=['Name','Highway','Amenity','geometry'],crs={'init': 'epsg:4326'})

    data_list = []
    for feature in lines:
        data=feature.ExportToJson(as_object=True)
        points = data['geometry']['coordinates']
        name = data['properties']['name']
        data_list.append([name, points])
    mapa['ruas'] = gpd.GeoDataFrame(data_list,columns=['Name','geometry'],crs={'init': 'epsg:4326'})

    return mapa

def _json_object_hook(d): return namedtuple('X', d.keys())(*d.values())
def json2obj(data): return json.loads(data, object_hook=_json_object_hook)

def equationLine():
    points = [[-52.616555,-27.0893808],[-52.6165129, -27.0893678]]
    x_coords, y_coords = zip(*points)
    A = vstack([x_coords,ones(len(x_coords))]).T
    m, c = lstsq(A, y_coords)[0]
    y = m*-52.616666+c

def linhas(square, ponto):
    linhas = []
    for i, row in mapa['ruas'].iterrows():
        for coords in row['geometry']:
            temp = [coords[0], coords[1]]
            ponto = Ponto(temp)
            if((ponto.y > square['minLatitude'] and  ponto.y < square['maxLatitude']) and (ponto.x > square['minLongitude'] and ponto.x < square['maxLongitude'])):
                linhas.append(row)
                break
    return linhas

def getData():
    file = sys.argv[1]
    #Open data file
    data=open(file).read()
    mapa = json2obj(data)
    mapa = json2obj(mapa.data)
    return mapa.locations

def plotMap(mapa):
    fig, ax = plt.subplots(figsize=(10,10))
    for row in mapa['pontos']['geometry']:
        x=row.x
        y=row.y
        plt.annotate('-', xy=(x,y), size=0.2, xytext=(0,5))
        plt.plot(x,y, '.', color='g')
        ax.set(aspect=1)

def plotTracos(tracos):
    for row in tracos:
        x = row.coords.longitude
        y = row.coords.latitude
        plt.annotate('-', xy=(x,y), size=2, xytext=(0,5))
        plt.plot(x,y, '.', color='r')

def plotSaida(dados):
    for row in dados:
        x = row[0]
        y = row[1]
        plt.annotate('-', xy=(x,y), size=2, xytext=(0,5))
        plt.plot(x,y, 'x', color= 'b')


def linkSelection(row, ponto):
    square = mm.square(row.coords.latitude, row.coords.longitude)
    candidatos = linhas(square, Ponto(ponto))
    pesos = []
    for candidato in candidatos:
        candidatosSorted = mm.getPointsNear(ponto, candidato['geometry'])
        #Pontos no link, mais próxmo ao ponto do gps
        A = candidatosSorted[0][0]
        saidaNear.append(A)
        #Ponto no link, segundo mais próximo do gps
        B = candidatosSorted[1][0]
        if(candidatosSorted[1][0] == candidatosSorted[0][0]):
            B = candidatosSorted[2][0]
        #Ponto mais ao norte
        C = [A[0],(A[1] +1)]
        deltaLinha = mm.deltaB(row.coords.heading, mm.getAngle(A, C, B, A))
        D = mm.D(Ponto(A), Ponto(B), Ponto(ponto))
        alpha = mm.getAngle(A, B, ponto, A)
        WSh = mm.WSh(np.cos(deltaLinha))
        WSpd = mm.WSpd(D)
        WSpi = mm.WSpi(mm.getAngle(A, B, ponto, oldPoint))
        WSrp = mm.WSrp(alpha)
        tws = mm.TWS(WSh, WSpd, WSpi, WSrp)
        pesos.append([candidato, tws, alpha])

    pesos = sorted(pesos, key=lambda item: item[1])
    return pesos[-1]

#Main
carregamentoInicio = time.time()
mapa = getMap()
data = getData()
carregamentoFim = time.time()
heading = -1
saidaMM = []
saidaMMNear = []
saidaNear = []
mm = GMM()
c=0
pE = 0
pN = 0
inicio = time.time()
selected = []

for row in data:
    ponto = [row.coords.longitude, row.coords.latitude]
    heading = row.coords.heading
    alpha = 0
    if(c>0):
        selectedLink = linkSelection(row, ponto) 
        selected = selectedLink[0]
        alpha = selectedLink[2]
        # Posição do veículo no segmento
        #ponto mais ao norte
        '''
        pn = [oldPoint[0], oldPoint[1]+1]
        theta = mm.getAngle(oldPoint, ponto, pn, oldPoint)
        deltaN = mm.deltaN(row.coords.speed, theta)
        deltaE = mm.deltaE(row.coords.speed, theta)
        pE = pE + deltaE
        pN = pN + deltaN
        possiblePoint = Ponto([pE, pN])
        '''
        possiblePoint = selected["geometry"][0]
        saidaMM.append(possiblePoint)
    else:
        pE = ponto[0]
        pN = ponto[1]

    oldPoint = ponto
    oldHeading = heading
    oldAlpha = alpha
    c+=1

fim = time.time()
plotMap(mapa)
plotTracos(data)
plotSaida(saidaMM)
total = fim - inicio
print("tempo(segundos): "+str(total))
print("carregamento mapas e coordenadas(segundos): "+ str(carregamentoFim - carregamentoInicio))
plt.show()