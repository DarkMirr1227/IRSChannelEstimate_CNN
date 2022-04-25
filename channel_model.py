from cmath import sqrt
from cmath import exp
from ctypes import sizeof
from random import random 
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d #3d 그래프 만들 때 사용
global wavelength,antenaGain,elementSizeM,elementSizeN,expone
#position값은 3차원 array라고 가정, 입사각인 theta도 구할 수 있게 프로그래밍해야함 
def positionSet(BSpos,IRSpos,UEpos):
    ''' bs,irs,ue의 3차원 포지션값을 이용하여 거리를 구하는 함수.
    
        args:
            BSpos: basestation의 (x,y,z) 값
            IRSpos: IRS의 정중앙(x,y,z) 값
            UEpos: UE의 (x,y,z) 값
        
        returns:
            dict('D_mn':bs과irs의 거리, 'd_mn':irs와ue의 거리, 'd_bu':bs와ue의 직선거리)
    '''
    D_mn = 0.
    d_mn = 0.
    d_bu = 0.
    #BS-IRS distance(D_mn) , IRS-UE distane(d_mn) , BS-UE direct distance(d_bu)
    for i in range(0,3,1):
        D_mn += (BSpos[i]-IRSpos[i])**2 
        d_mn += (IRSpos[i]-UEpos[i])**2
        d_bu += (BSpos[i]-UEpos[i])**2
    D_mn = math.sqrt(D_mn)
    d_mn = math.sqrt(d_mn)
    d_bu = math.sqrt(d_bu)
    return {'D_mn':D_mn,'d_mn':d_mn,'d_bu':d_bu}

def calcTheta(BSpos,IRSpos,UEpos):
    '''bs,irs,ue의 3차원 포지션값을 이용하여 입사각을 구하는 함수.

        args:
            BSpos: basestation의 (x,y,z) 값
            IRSpos: IRS의 정중앙(x,y,z) 값
            UEpos: UE의 (x,y,z) 값

        return:
            입사각(theta) : 세점 사이의 각도의 절반        
    '''
    v1 = BSpos-IRSpos
    v2 = UEpos-IRSpos
    return np.inner(v1,v2)/(np.linalg.norm(v1) *np.linalg.norm(v2)*2)

def initValue(_wavelength,_antenaGain=1,_elementSizeM=0.04,_elementSizeN=0.04,_expone=2):
    ''' 수식계산을 위한 초기값(상수)을 설정하는 함수.
        
    입력값을 글로벌변수에 대입한다.
    
        args:
            _wavelength: 파장길이( m단위)
            _antenaGain: 안테나 이득(default=1)
            _elementSizeM: IRS element의 세로 길이
            _elementSizeN: IRS element의 가로 길이
            _expone: 신호를 구하는 수식에서의 alpha값 (default=2)
        
        return:
            없음
    '''
    global wavelength,antenaGain,elementSizeM,elementSizeN,expone
    wavelength = _wavelength
    antenaGain = _antenaGain
    elementSizeM = _elementSizeM
    elementSizeN = _elementSizeM
    expone = _expone

def calc(_theta,_D_mn,_d_mn,_d_bu):
    '''변수값을 입력받아 r_mn*h_mn을 구하는 함수

        args:
            _theta: 입사각
            _D_mn: BS와 IRS의 거리
            _d_mn: IRS와 UE의 거리
            _d_bu: BS와 UE의 직선거리
        
        return:
            r_mn*h_mn
    '''
    global wavelength,antenaGain,elementSizeM,elementSizeN,expone
    result1 = math.cos(_theta) * wavelength *sqrt(antenaGain*elementSizeM*elementSizeN) / ((4*math.pi)**(3/2)*sqrt(_D_mn*_d_mn)**expone)
    temp = (math.pi * _d_bu / wavelength)%(2*math.pi)
    result2 = exp(-2j * temp)
    return result1*result2

def setElementPosXZ(IRSpos,mNum,nNum,delta=0):
    '''IRS 각각의 element position값을 구하는 함수. 여기서 IRS의 element의 Y값은 모두 동일한 경우에만 사용가능.

    즉, 변하는값이 "x","z"일때만 해당 함수를 사용할 수 있음.

        args:
            IRSpos: IRS의 정중앙(x,y,z) 값
            mNum: IRS의 총 element 가로 갯수
            nNum: IRS의 총 element 세로 갯수
            delta: element 사이의 간격/distance (default =0)
        
        return:
            elementPos: 각 element(m,n)에 따른 position값
    '''
    elementPos = np.zeros((mNum,nNum,3))
    for m in range(0,mNum,1):
        for n in range(0,nNum,1):
            elementPos[m][n][0] = IRSpos[0] -(mNum+1)/2*(delta+elementSizeN) + (n+1)*(delta+elementSizeN)
            elementPos[m][n][1] = IRSpos[1]
            elementPos[m][n][2] = IRSpos[2] +(nNum+1)/2*(delta+elementSizeM) - (m+1)*(delta+elementSizeM)

    return elementPos

def allElementCalc(BSpos,IRSpos,UEpos,mNum,nNum):
    '''각 IRS element 수만큼 거리,theta값을 구하여 r_mn*h_mn을 구하는 함수
    
        args:
            BSpos: basestation의 (x,y,z) 값
            IRSpos: IRS의 정중앙(x,y,z) 값
            UEpos: UE의 (x,y,z) 값
            mNum: IRS의 총 element 가로 갯수
            nNum: IRS의 총 element 세로 갯수

        return:
            elementResult: 각 element(m,n)의 r_mn*h_mn 값            
            
    '''
    elementPos = setElementPosXZ(IRSpos,mNum,nNum)
    elementResult = np.zeros((mNum,nNum),dtype=np.complex64)
    for m in range(0,mNum,1):
        for n in range(0,nNum,1):
            dicDist = positionSet(BSpos,elementPos[m][n],UEpos)
            theta = calcTheta(BSpos,elementPos[m][n],UEpos)
            elementResult[m][n] = calc(theta,dicDist['D_mn'],dicDist['d_mn'],dicDist['d_bu'])

    return elementResult
 
initValue(0.01)
Baseposition = np.array([0,50,35]) # 0,50,35
IRSposition = np.array([75,100,2])
UEposition = np.array([50+random()*50,random()*50,1.+random()])

result = allElementCalc(Baseposition,IRSposition,UEposition,mNum=6,nNum=6)
# dicDist = positionSet(Baseposition,IRSposition,UEposition)
# theta = calcTheta(Baseposition,IRSposition,UEposition)
# temp = calc(theta,dicDist['D_mn'],dicDist['d_mn'],dicDist['d_bu'])

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
X = np.linspace(0,5,6,dtype=np.int64)
Y = np.linspace(0,5,6,dtype=np.int64)
X,Y = np.meshgrid(X,Y)
ax.plot_surface(X,Y,result)
#print(result)
print('UEposition : ',UEposition)
plt.show()

#datalist에 data 쌓기
datalist = []
for count in range(0,100,1):
    UEposition = np.array([50+random()*50,random()*50,1.5+random()-0.5])
    result = allElementCalc(Baseposition,IRSposition,UEposition,mNum=6,nNum=6)
    datalist.append(result)

print(len(datalist))

#CNN

