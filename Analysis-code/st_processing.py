import math
import numpy as np

def make_uniq_list(a):
    a=list(a)
    return [ str(a[i]) + '_' + str(a[0:i].count(a[i])+1) if a[0:i].count(a[i]) > 0 else a[i] for i in range(len(a))]

def load_image(filename):
    print(f'Image loading from {filename}')
    from PIL import Image
    from PIL import ImageFile
    ImageFile.LOAD_TRUNCATED_IMAGES = True
    Image.MAX_IMAGE_PIXELS = None
    img = Image.open(filename)
    img = np.array(img)
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
        # remove alpha channel
    return img
def rgba_to_hex(rgba):
    r, g, b, a = rgba
    # 将RGB分量从0到1的范围转换为0到255的范围
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    a =int(a * 100)
    if a==100:
        a=''
    # 将RGB分量转换为HEX格式
    hex_color = "#{:02X}{:02X}{:02X}{}".format(r, g, b,a)
    return hex_color
def hex2rgb(hex_str):
    h = hex_str.lstrip('#')
    rgb = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    r,g,b=rgb
    rgb=(r/255,g/255,b/255)
    return rgb

#return grey image
def rgb2gray(rgb):
    import numpy as np
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140]).astype('uint8') 

# 根据已知两点坐标，求过这两点的直线解析方程： a*x+b*y+c = 0  (a >= 0)
def getLinearEquation_new(line_point1,line_point2):
    p1x=line_point1[0]
    p1y=line_point1[1]
    p2x=line_point2[0]
    p2y=line_point2[1]
    k=(p1y-p2y)/(p1x-p2x)*1.0
    # if k < 0:
    #     k = -k
    b = ((p1y+p2y) - k * (p1x+p2x))/2
    return [k, -1, b]

#一个点到一条直线的距离
def get_distance_from_point_to_line(point, a, b, c):
    import numpy as np
    #a*x+b*y+c = 0  (a >= 0)
    #计算直线的三个参数
    A = a
    B = b
    C = c
    #根据点到直线的距离公式计算距离
    distance = np.abs(A * point[0] + B * point[1] + C) / (np.sqrt(A**2 + B**2))
    return [A, B, C ,distance]

def dis_two_point(point1,point2):
    import math
    x=point1[0]-point2[0]
    y=point1[1]-point2[1]
    #用math.sqrt（）求平方根
    return(math.sqrt((x**2)+(y**2)))

def get_piexl_coods(line_point_r1c1,line_point_r70c1,line_point_r1c70,barcode_num,reg,chip_type):
    import numpy as np
    import math
    x1_a,x1_b,x1_c=getLinearEquation_new(line_point_r1c1,line_point_r70c1)#col line_A1
    y1_a,y1_b,y1_c=getLinearEquation_new(line_point_r1c1,line_point_r1c70)#row line_B1

    #第一个点到line2的距离
    x_dis=get_distance_from_point_to_line(line_point_r1c70, x1_a,x1_b,x1_c)[3]/(barcode_num-1)

    y_dis=get_distance_from_point_to_line(line_point_r70c1,y1_a,y1_b,y1_c)[3]/(barcode_num-1)

    line_X_A={}
    line_Y_B={}

    n=0
    for i in range(0,barcode_num):
        n+=1
        if x1_a < 0 :
            c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*i,2) + pow(x_dis*i*(x1_a/x1_b),2)) 
            line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
        else:
            c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*i,2) + pow(x_dis*i*(x1_a/x1_b),2)) 
            line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
        if y1_a < 0 :
            c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*i,2) + pow(y_dis*i*(y1_a/y1_b),2)) 
            line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
        else:
            c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*i,2) + pow(y_dis*i*(y1_a/y1_b),2)) 
            line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)
        i=i+1
    pixel_coords=[]
    
    if (chip_type=='T9') and (reg in ['reg1','reg3','reg7','reg9']):
        for m in range(barcode_num,0,-1):
            for n in range(1,(barcode_num+1)):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords.append([int(x_col),int(y_row)])
    elif (chip_type=='T9') and (reg in ['reg2','reg8']):
        for m in range(1,(barcode_num+1)):
            for n in range(1,(barcode_num+1)):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords.append([int(x_col),int(y_row)])
    elif (chip_type=='T9') and (reg in ['reg4','reg6']):
        for m in range(barcode_num,0,-1):
            for n in range(barcode_num,0,-1):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords.append([int(x_col),int(y_row)])  
    elif (chip_type=='T9') and (reg in ['reg5']):
        for m in range(1,(barcode_num+1)):
            for n in range(barcode_num,0,-1):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords.append([int(x_col),int(y_row)]) 
    elif (chip_type=='T3') and (reg in ['reg1','reg3']):
        for m in range(1,(barcode_num+1)):
            for n in range(1,(barcode_num+1)):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords.append([int(x_col),int(y_row)])
    elif (chip_type=='T3') and (reg in ['reg2']):
        for m in range(barcode_num,0,-1):
            for n in range(1,(barcode_num+1)):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords.append([int(x_col),int(y_row)])
    else:
        print('No chip coord')
    return np.array(pixel_coords),(x_dis+y_dis)/4
def get_HE_mask(img_HE):
    import numpy as np
    from skimage import data,filters
    import cv2
    image=rgb2gray(img_HE)
    thresh = filters.threshold_otsu(image)   #返回一个阈值

    _, image_b = cv2.threshold(image, thresh, 250, cv2.THRESH_BINARY_INV)
    mask=np.zeros([image.shape[0],image.shape[1]])
    # 找到所有的轮廓
    contours, _ = cv2.findContours(image_b, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    area = []
    # 找到最大的轮廓
    for k in range(len(contours)):
        area.append(cv2.contourArea(contours[k]))
    max_idx = np.argmax(np.array(area))
    # 填充最大的轮廓
    mask = cv2.drawContours(mask, contours, max_idx, 1, cv2.FILLED)
    return mask

def get_HE_mask_3d(img_file,img,image_file_path):
    import numpy as np
    from skimage import data,filters
    import cv2
    dx=2.98
    if img_file.startswith(image_file_path+'T9-5-') or img_file.startswith(image_file_path+'T9-11-'):
        img = cv2.convertScaleAbs(img, alpha=1.5, beta=0)#增强对比度
    if img_file.startswith(image_file_path+'T9-5-') or img_file.startswith(image_file_path+'T9-11-') or img_file.startswith(image_file_path+'T9-17-7'):
        img[0:int(2600/dx),:]=255
        img[:,0:int(2300/dx)]=255
        img[:,int(12000/dx):]=255
        img[int(11000/dx):,:]=255
    elif  img_file.startswith(image_file_path+'T9-8-7'):
        img[0:int(1600/dx),:]=255
        img[:,0:int(2300/dx)]=255
        img[:,int(12000/dx):]=255
        img[int(11000/dx):,:]=255
    elif  img_file.startswith(image_file_path+'T9-10-6'):
        img[0:int(1600/dx),:]=255
        img[:,0:int(2300/dx)]=255
        img[:,int(12000/dx):]=255
        img[int(13300/dx):,:]=255
    else:    
        img[0:int(2600/dx),:]=255
        img[:,0:int(2300/dx)]=255
        img[:,int(12000/dx):]=255
        img[int(12800/dx):,:]=255
    if len(img.shape)>2:
        image=rgb2gray(img)
    else:
        image=img
    thresh = filters.threshold_otsu(image)   #返回一个阈值

    _, image_b = cv2.threshold(image, thresh, 250, cv2.THRESH_BINARY_INV)
    mask=np.zeros([image.shape[0],image.shape[1]])
    mask_tmp=np.zeros([image.shape[0],image.shape[1]])
    # 找到所有的轮廓
    contours, _ = cv2.findContours(image_b, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    area = []
    # 找到最大的轮廓
    for k in range(len(contours)):
        area.append(cv2.contourArea(contours[k]))
    max_idx = np.argmax(np.array(area))
    mask = cv2.drawContours(mask, contours, max_idx, 1, cv2.FILLED)
    if np.array(area).shape[0]>1:
        idx = np.argsort(np.array(area))[-2]
    # 填充最大的轮廓
        mask_tmp = cv2.drawContours(mask_tmp, contours, idx, 1, cv2.FILLED)
        if mask_tmp.sum() > mask.sum()/10:
            mask = cv2.drawContours(mask, contours, idx, 1, cv2.FILLED)
    return mask
def get_chip_coord(barcode_num,reg,chip_type,Barcode_file_path=None):
    import pandas as pd
    file_A1_70=Barcode_file_path+chip_type+'-A1-A'+str(barcode_num)+'.txt'
    file_B1_70=Barcode_file_path+chip_type+'-B1-B'+str(barcode_num)+'.txt'

    A1_70=['A' + str(i) for i in pd.read_csv(file_A1_70,header=None)[0]]
    B1_70=['B' + str(i) for i in pd.read_csv(file_B1_70,header=None)[0]]

    if (chip_type=='T9') and (reg in ['reg1','reg3','reg7','reg9']):
        df_A1_70=pd.DataFrame({'col':list(range(barcode_num,0,-1))},index=A1_70)
        df_B1_70=pd.DataFrame({'row':list(range(1,barcode_num+1,1))},index=B1_70)
    elif (chip_type=='T9') and (reg in ['reg2','reg8']):
        df_A1_70=pd.DataFrame({'col':list(range(1,barcode_num+1,1))},index=A1_70)
        df_B1_70=pd.DataFrame({'row':list(range(1,barcode_num+1,1))},index=B1_70)
    elif (chip_type=='T9') and (reg in ['reg4','reg6']):
        df_A1_70=pd.DataFrame({'col':list(range(barcode_num,0,-1))},index=A1_70)
        df_B1_70=pd.DataFrame({'row':list(range(barcode_num,0,-1))},index=B1_70)    
    elif (chip_type=='T9') and (reg in ['reg5']):
        df_A1_70=pd.DataFrame({'col':list(range(1,barcode_num+1,1))},index=A1_70)
        df_B1_70=pd.DataFrame({'row':list(range(barcode_num,0,-1))},index=B1_70)  
    elif (chip_type=='T3') and (reg in ['reg2']):
        df_A1_70=pd.DataFrame({'col':list(range(1,barcode_num+1,1))},index=A1_70)
        df_B1_70=pd.DataFrame({'row':list(range(1,barcode_num+1,1))},index=B1_70) 
    elif (chip_type=='T3') and (reg in ['reg1','reg3']):
        df_A1_70=pd.DataFrame({'col':list(range(barcode_num,0,-1))},index=A1_70)
        df_B1_70=pd.DataFrame({'row':list(range(1,barcode_num+1,1))},index=B1_70) 
    else:
        print('No chip coord')

    Spot_coordinate={'col':[],
    'row':[]}
    Spot=[]
    for a in df_A1_70.index:
        for b in df_B1_70.index:
            Spot.append(a+'-'+b)
            Spot_coordinate['col'].append(df_A1_70.loc[a,'col'])
            Spot_coordinate['row'].append(df_B1_70.loc[b,'row'])
    Spot_coordinate=pd.DataFrame.from_dict(Spot_coordinate,orient='index',columns=Spot).T
    return Spot_coordinate
def get_piexl_coods_Mchip(line_point_r1c1,line_point_r210c1,line_point_r1c210,channels_num,barcode_num,chip_type,Barcode_file_path,add_M9_15=False,add_M9_20=False):
    file_C1_9=Barcode_file_path+chip_type+'-C1-C'+str(18)+'.txt'
    x1_a,x1_b,x1_c=getLinearEquation_new(line_point_r1c1,line_point_r210c1)#col line_A1
    y1_a,y1_b,y1_c=getLinearEquation_new(line_point_r1c1,line_point_r1c210)#row line_B1
    import numpy as np
    import pandas as pd
    #第一个点到line2的距离
    if barcode_num==70 and chip_type=='M9':
        x_dis=get_distance_from_point_to_line(line_point_r1c210, x1_a,x1_b,x1_c)[3]/((barcode_num-1)*3+7)
        y_dis=get_distance_from_point_to_line(line_point_r210c1,y1_a,y1_b,y1_c)[3]/((barcode_num-1)*3+7)
    elif barcode_num==150 and chip_type=='M9':
        x_dis=get_distance_from_point_to_line(line_point_r1c210, x1_a,x1_b,x1_c)[3]/((barcode_num-1)*3+15)
        y_dis=get_distance_from_point_to_line(line_point_r210c1,y1_a,y1_b,y1_c)[3]/((barcode_num-1)*3+15)
    else:
        print(f'Chip type error: {chip_type}')
    line_X_A={}
    line_Y_B={}

    n=0
    if barcode_num==70 and chip_type=='M9':
        for i in range(0,channels_num):
            n+=1
            if n>0 and n<(barcode_num+1) :
                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*i,2) + pow(x_dis*i*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*i,2) + pow(x_dis*i*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*i,2) + pow(y_dis*i*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*i,2) + pow(y_dis*i*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)

            elif n>barcode_num and n<(barcode_num*2+1) :
                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+2.5),2) + pow(x_dis*(i+2.5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+2.5),2) + pow(x_dis*(i+2.5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*(i+2.5),2) + pow(y_dis*(i+2.5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*(i+2.5),2) + pow(y_dis*(i+2.5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)
            elif n>barcode_num*2 and n<(barcode_num*3+1) :
                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+5),2) + pow(x_dis*(i+5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+5),2) + pow(x_dis*(i+5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*(i+5),2) + pow(y_dis*(i+5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*(i+5),2) + pow(y_dis*(i+5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)    
    elif barcode_num==150 and chip_type=='M9':
        for i in range(0,channels_num):
            n+=1
            if n>0 and n<(barcode_num+1) :
                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*i,2) + pow(x_dis*i*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*i,2) + pow(x_dis*i*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*i,2) + pow(y_dis*i*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*i,2) + pow(y_dis*i*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)

            elif n>barcode_num and n<(barcode_num*2+1) :
                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+6.5),2) + pow(x_dis*(i+6.5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+6.5),2) + pow(x_dis*(i+6.5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*(i+6.5),2) + pow(y_dis*(i+6.5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*(i+6.5),2) + pow(y_dis*(i+6.5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)
            elif n>barcode_num*2 and n<(barcode_num*3+1) :
                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+13),2) + pow(x_dis*(i+13)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+13),2) + pow(x_dis*(i+13)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*(i+13),2) + pow(y_dis*(i+13)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*(i+13),2) + pow(y_dis*(i+13)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)    
    else:
        print(f'Chip type error: {chip_type}')
    if add_M9_15:
        for i in [barcode_num ,barcode_num+1 ,barcode_num*2 ,barcode_num*2+1]:
            n+=1
            if i < (barcode_num+2):

                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+0.25),2) + pow(x_dis*(i+0.25)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+0.25),2) + pow(x_dis*(i+0.25)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*(i+0.25),2) + pow(y_dis*(i+0.25)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*(i+0.25),2) + pow(y_dis*(i+0.25)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)
            elif i > (barcode_num+2) and i < (barcode_num*2+2) :

                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+0.25+2.5),2) + pow(x_dis*(i+0.25+2.5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+0.25+2.5),2) + pow(x_dis*(i+0.25+2.5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*(i+0.25+2.5),2) + pow(y_dis*(i+0.25+2.5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*(i+0.25+2.5),2) + pow(y_dis*(i+0.25+2.5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2) 
        all_n=n
    elif add_M9_20:
        for i in [barcode_num ,barcode_num+1 ,barcode_num+2 ,barcode_num+3 ,barcode_num+4 ,barcode_num+5 ,barcode_num+6,
        barcode_num*2 ,barcode_num*2+1,barcode_num*2+2,barcode_num*2+3,barcode_num*2+4,barcode_num*2+5,barcode_num*2+6]:
            n+=1
            if i < (barcode_num+7):

                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*i,2) + pow(x_dis*i*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*i,2) + pow(x_dis*i*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*i,2) + pow(y_dis*i*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*i,2) + pow(y_dis*i*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)
            elif i > (barcode_num+7) and i < (barcode_num*2+7) :

                if x1_a < 0 :
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+6.5),2) + pow(x_dis*(i+6.5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)     
                else:
                    c1=x1_c+abs(x1_b)*math.sqrt(pow(x_dis*(i+6.5),2) + pow(x_dis*(i+6.5)*(x1_a/x1_b),2)) 
                    line_X_A['X_A_line'+str(n)]=(x1_a,x1_b,c1)
                if y1_a < 0 :
                    c2=y1_c+abs(y1_b)*math.sqrt(pow(y_dis*(i+6.5),2) + pow(y_dis*(i+6.5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2)     
                else:
                    c2=y1_c-abs(y1_b)*math.sqrt(pow(y_dis*(i+6.5),2) + pow(y_dis*(i+6.5)*(y1_a/y1_b),2)) 
                    line_Y_B['Y_B_line'+str(n)]=(y1_a,y1_b,c2) 
        all_n=n

    pixel_coords=[]
    C1_9=['C' + str(i) for i in pd.read_csv(file_C1_9,header=None)[0]][0:18:2]
    #C1_9=list(Spot_coordinate_uniq.drop_duplicates(['barcode-C'])['barcode-C'])
    for reg_i in C1_9:
        if (chip_type=='M9') and (reg_i == C1_9[0]):
            for m in range(barcode_num,0,-1):
                for n in range(1,(barcode_num+1)):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)])
        elif (chip_type=='M9') and (reg_i == C1_9[1]):
            for m in range(1+barcode_num,(barcode_num+1+barcode_num)):
                for n in range(1,(barcode_num+1)):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)])
        elif (chip_type=='M9') and (reg_i == C1_9[2]):
            for m in range(barcode_num+barcode_num*2,0+barcode_num*2,-1):
                for n in range(1,(barcode_num+1)):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)])                
        elif (chip_type=='M9') and (reg_i == C1_9[3]):
            for m in range(barcode_num,0,-1):
                for n in range(barcode_num+barcode_num,0+barcode_num,-1):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)])  
        elif (chip_type=='M9') and (reg_i == C1_9[4]):
            for m in range(1+barcode_num,(barcode_num+1+barcode_num)):
                for n in range(barcode_num+barcode_num,0+barcode_num,-1):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)]) 
        elif (chip_type=='M9') and (reg_i == C1_9[5]):
            for m in range(barcode_num+barcode_num*2,0+barcode_num*2,-1):
                for n in range(barcode_num+barcode_num,0+barcode_num,-1):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)])    
        elif (chip_type=='M9') and (reg_i == C1_9[6]):
            for m in range(barcode_num,0,-1):
                for n in range(1+barcode_num*2,(barcode_num+1+barcode_num*2)):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)])    
        elif (chip_type=='M9') and (reg_i == C1_9[7]):
            for m in range(1+barcode_num,(barcode_num+1+barcode_num)):
                for n in range(1+barcode_num*2,(barcode_num+1+barcode_num*2)):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)])
        elif (chip_type=='M9') and (reg_i == C1_9[8]):
            for m in range(barcode_num+barcode_num*2,0+barcode_num*2,-1):
                for n in range(1+barcode_num*2,(barcode_num+1+barcode_num*2)):
                    a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                    a1,b1,c1=line_X_A['X_A_line'+str(m)]
                    D = a0*b1 - a1*b0
                    y_row = (b0*c1 - b1*c0)/D
                    x_col = (a1*c0 - a0*c1)/D
                    pixel_coords.append([int(x_col),int(y_row)])
        else:
            print('No chip coord')
            
    if add_M9_15:
        pixel_coords_all=[]
        spot_coords_all=[]
        for m in [channels_num + i for i in range(1,5)]:
            for n in range(1,(all_n+1)):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords_all.append([int(x_col),int(y_row)])
                spot_coords_all.append([int(m),int(n)])
        for n in [channels_num + i for i in range(1,5)]:
            for m in range(1,(all_n+1)):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords_all.append([int(x_col),int(y_row)])
                spot_coords_all.append([int(m),int(n)])
        pixel_coords_all=np.array(list_remove_dup(pixel_coords_all)).astype(int)
        spot_coords_all=np.array(list_remove_dup(spot_coords_all)).astype(int)
        return np.array(pixel_coords),pixel_coords_all,spot_coords_all,(x_dis+y_dis)/4
    elif add_M9_20:
        pixel_coords_all=[]
        spot_coords_all=[]
        for m in [channels_num + i for i in range(1,15)]:
            for n in range(1,(all_n+1)):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords_all.append([int(x_col),int(y_row)])
                spot_coords_all.append([int(m),int(n)])
        for n in [channels_num + i for i in range(1,15)]:
            for m in range(1,(all_n+1)):
                a0,b0,c0=line_Y_B['Y_B_line'+str(n)]
                a1,b1,c1=line_X_A['X_A_line'+str(m)]
                D = a0*b1 - a1*b0
                y_row = (b0*c1 - b1*c0)/D
                x_col = (a1*c0 - a0*c1)/D
                pixel_coords_all.append([int(x_col),int(y_row)])
                spot_coords_all.append([int(m),int(n)])
        pixel_coords_all=np.array(list_remove_dup(pixel_coords_all)).astype(int)
        spot_coords_all=np.array(list_remove_dup(spot_coords_all)).astype(int)
        return np.array(pixel_coords),pixel_coords_all,spot_coords_all,(x_dis+y_dis)/4
    else:
        return np.array(pixel_coords),(x_dis+y_dis)/4  # x_dis 两个spot中心点距离


    
def list_remove_dup(input_list):
    new_list = []
    for element in input_list:
        if element not in new_list:
            new_list.append(element)
    return new_list

def trans_pN(x):
    x=np.array(x,dtype=np.float32)
    P_99=np.nanpercentile(x, 99)
    P_1=np.nanpercentile(x, 1)
    x[x>P_99]=P_99
    x[x<P_1]=P_1
    return x 
    
###计算各种分位数
def calculateStuff4Sam(i):
    import numpy as np
    res = {}
    numbers = np.array(i)
    res['mean'] = "{0:.2f}".format(np.mean(numbers))
    res['median'] = "{0:.2f}".format(np.median(numbers))
    res['min'] = "{0:.2f}".format(np.min(numbers))
    res['25 percentile'] = "{0:.2f}".format(np.percentile(numbers, 25))
    res['50 percentile'] = "{0:.2f}".format(np.percentile(numbers, 50))
    res['75 percentile'] = "{0:.2f}".format(np.percentile(numbers, 75))
    res['max'] = "{0:.2f}".format(np.max(numbers))
    return res
#打印重复的值
def return_dup(input_list):
    lst = list(input_list)
    duplicates = set()
    for i in lst:
        if lst.count(i) > 1:
            duplicates.add(i)
    print(duplicates)

    
def get_hex_from_cmap(tab10):
    colors = []
    for i in range(tab10.shape[0]):
        r=int(256*tab10[i][0])
        g=int(256*tab10[i][1])
        b=int(256*tab10[i][2])
        hex_color = f'#{r:02x}{g:02x}{b:02x}'
        colors.append(hex_color)
    return colors
#tab10=cm.get_cmap('tab10',10).colors[:,0:3]
#get_hex_from_cmap(tab10)

#按指定点旋转坐标点-逆时针旋转
def Nrotate(angle, point, centre_point):
    x,y=point
    cx,cy=centre_point
    angle = math.radians(angle)
    x = np.array(x)
    y = np.array(y)
    nRotatex = (x - cx) * math.cos(angle) - (y - cy) * math.sin(angle) + cx
    nRotatey = (x - cx) * math.sin(angle) + (y - cy) * math.cos(angle) + cy
    return round(nRotatex, 2), round(nRotatey, 2)
#按指定点旋转坐标点-顺时针旋转 
def Srotate(angle, point, centre_point):
    x,y=point
    cx,cy=centre_point
    angle = math.radians(angle)
    x = np.array(x)
    y = np.array(y)
    sRotatex = (x - cx) * math.cos(angle) + (y - cy) * math.sin(angle) + cx
    sRotatey = (y - cy) * math.cos(angle) - (x - cx) * math.sin(angle) + cy
    return round(sRotatex, 2), round(sRotatey, 2)

#Adata-顺时针旋转 
def adata_Srotate(adata,angle=0,borderValue=(0,0,0)):
    import cv2 
    import numpy as np
    adata=adata.copy()
    spatial_coord=adata.obsm['spatial'].copy()#（col,row）
    spatial_coord=spatial_coord[:,::-1]#（row,col）
    sample_id=list(adata.uns['spatial'].keys())[0]
    img_hires=adata.uns['spatial'][sample_id]['images']['hires']
    img_lowres=adata.uns['spatial'][sample_id]['images']['lowres']
    ths=adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']
    tls=adata.uns['spatial'][sample_id]['scalefactors']['tissue_lowres_scalef']
    n_col=img_hires.shape[1]/ths
    n_row=img_hires.shape[0]/ths
    centre_point=(n_row/2,n_col/2)#(rows,cols)
    #按指定点旋转坐标点-顺时针旋转 
    for r in range(spatial_coord.shape[0]):
        point=spatial_coord[r,:]
        point_new=Srotate(angle, point, centre_point)
        spatial_coord[r,:]=point_new
    spatial_coord=spatial_coord.astype(int)
    
    if int(np.min(spatial_coord))<0:
        m_m=abs(int(np.min(spatial_coord)))
        spatial_coord=spatial_coord+m_m
    else:
        m_m=0

    adata.obsm['spatial']=spatial_coord[:,::-1]#（col,row）
    ##img_hires
    #旋转
    m=int(np.max(spatial_coord)*ths)#max(img_hires.shape[0],img_hires.shape[1])
    centre_point=(ths*n_col/2,ths*n_row/2) # (cols,rows)
    M = cv2.getRotationMatrix2D(centre_point,-angle,1)
    img_hires = cv2.warpAffine(img_hires,M,(m,m),borderValue=borderValue)
    #平移
    mm=int(m_m*ths)
    M = np.float32([[1, 0, mm], [0, 1, mm]]) #第1个 左右平移  第2个 上下平移
    # 用仿射变换实现平移
    img_hires = cv2.warpAffine(img_hires, M, (m,m),borderValue=borderValue) #cols, row 正常的 row col

    #img_lowres
    m=int(np.max(spatial_coord)*tls)#max(img_lowres.shape[0],img_lowres.shape[1])
    centre_point=(tls*n_col/2,tls*n_row/2) # (cols,rows)
    M = cv2.getRotationMatrix2D(centre_point,-angle,1)
    img_lowres = cv2.warpAffine(img_lowres,M,(m,m),borderValue=borderValue)
    #平移
    mm=int(m_m*ths)
    M = np.float32([[1, 0, mm], [0, 1, mm]]) #第1个 左右平移  第2个 上下平移
    # 用仿射变换实现平移
    img_lowres = cv2.warpAffine(img_lowres, M, (m,m),borderValue=borderValue) #cols, row 正常的 row col
    
    adata.uns['spatial'][sample_id]['images']['hires']=img_hires
    adata.uns['spatial'][sample_id]['images']['lowres']=img_lowres
    return adata
    
#Adata-翻转
def adata_flip(adata,flip=0):#flip 0 --- 垂直方向翻转； 1----- 水平方向翻转； -1：水平、垂直方向同时翻转
    import cv2 
    adata=adata.copy()
    sample_id=list(adata.uns['spatial'].keys())[0]
    spatial_coord=adata.obsm['spatial'].copy()#（col,row）
    img_hires=adata.uns['spatial'][sample_id]['images']['hires']
    img_lowres=adata.uns['spatial'][sample_id]['images']['lowres']
    ths=adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']
    
    m_row=int(img_hires.shape[0]/ths)
    m_col=int(img_hires.shape[1]/ths)
    
        
    if flip==0:
        # Flipped Vertically 垂直翻转
        img_hires = cv2.flip(img_hires, 0)
        img_lowres = cv2.flip(img_lowres, 0)
        spatial_coord[:,1]=m_col-spatial_coord[:,1]
    elif flip==1:
        # Flipped Horizontally 水平翻转
        img_hires = cv2.flip(img_hires, 1)
        img_lowres = cv2.flip(img_lowres, 1)
        spatial_coord[:,0]=m_col-spatial_coord[:,0]
    elif flip==-1:
        # Flipped Horizontally & Vertically 水平垂直翻转
        img_hires = cv2.flip(img_hires, -1)
        img_lowres = cv2.flip(img_lowres, -1)
        spatial_coord[:,0]=m_col-spatial_coord[:,0]
        spatial_coord[:,1]=m_col-spatial_coord[:,1]
    else:
        print('No flip')
    
    adata.obsm['spatial']=spatial_coord
    adata.uns['spatial'][sample_id]['images']['hires']=img_hires
    adata.uns['spatial'][sample_id]['images']['lowres']=img_lowres
    return adata
    

    


def get_adata_STARsolo(sample,chip_type,reg,channels_num,barcode_num,res_um=None,only_HE=False,s_id=False,s=None,
            image_pre=False,# For 3D brain image preprocessing
            mask_sam=False,
            add_M9_15=False,
            add_M9_20=False,
            img_file=True,EM=True,Velocyto=False,soloFeatures='Gene',raw=True,species='human',
            file_path=None,image_file_path=None,Barcode_file_path=None,geneinfo_path=None,
            HE_point1=None,Spot_point1=None,HE_point2=None,Spot_point2=None,
            line_point_r1c1=None,line_point_r1c70=None,line_point_r70c1=None):
    import time
    start_time = time.time()
    import cv2 
    import scanpy as sc
    import os
    import numpy as np
    import pandas as pd
    from scipy.io import mmread
    from anndata import AnnData
    import scipy.sparse as sp
    print(f'sample: {sample}')
    if Velocyto:
        print('Reading Velocyto')         
        if raw:
            #x_ambiguous=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/raw/ambiguous.mtx.gz').astype(int).toarray().T)
            x_spliced=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/raw/spliced.mtx.gz').astype(int).toarray().T)
            x_unspliced=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/raw/unspliced.mtx.gz').astype(int).toarray().T)
        else:
            #x_ambiguous=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/filtered/ambiguous.mtx.gz').astype(int).toarray().T)
            x_spliced=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/filtered/spliced.mtx.gz').astype(int).toarray().T)
            x_unspliced=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/filtered/unspliced.mtx.gz').astype(int).toarray().T)
    else:
        print('No Velocyto')         
            
    if soloFeatures=='Gene' or soloFeatures=='GeneFull':
        if raw:
            if EM:
                gene_file=file_path+sample+'_Solo.out/'+soloFeatures+'/raw/UniqueAndMult-EM.mtx.gz'
            else:
                gene_file=file_path+sample+'_Solo.out/'+soloFeatures+'/raw/matrix.mtx.gz'
            barcodes_file=file_path+sample+'_Solo.out/'+soloFeatures+'/raw/barcodes.tsv'
            features_file=file_path+sample+'_Solo.out/'+soloFeatures+'/raw/features.tsv'
        else:
            gene_file=file_path+sample+'_Solo.out/'+soloFeatures+'/filtered/matrix.mtx.gz'
            barcodes_file=file_path+sample+'_Solo.out/'+soloFeatures+'/filtered/barcodes.tsv'
            features_file=file_path+sample+'_Solo.out/'+soloFeatures+'/filtered/features.tsv'
    else:
        print('No soloFeatures')

    geneinfo_path=geneinfo_path+'geneInfo.txt'

    if img_file:
        if s_id:
            sample_name=sample[:sample.rfind(s)]
        else:
            sample_name = sample
        HE_file=chip_type+'-'+sample_name+'-HE'
        HE_file=[i for i in os.listdir(image_file_path) if HE_file in i ]

        if len(HE_file)>0:
            HE_file=HE_file[0]
            HE_file=os.path.join(image_file_path,HE_file)
        else:
            raise FileNotFoundError('HE image not found')
            
        spot_marker_file=chip_type+'-'+sample_name+'-marker'
        spot_marker_file=[i for i in os.listdir(image_file_path) if spot_marker_file in i ]
        if len(spot_marker_file)>0:
            spot_marker_file=spot_marker_file[0]
            spot_marker_file=os.path.join(image_file_path,spot_marker_file)
        
        spot_file=chip_type+'-'+sample_name+'-spot'
        spot_file=[i for i in os.listdir(image_file_path) if spot_file in i ]
        if len(spot_file)>0:
            spot_file=spot_file[0]
            spot_file=os.path.join(image_file_path,spot_file)
        
        mask_file=chip_type+'-'+sample_name+'-mask'
        mask_file=[i for i in os.listdir(image_file_path) if mask_file in i ]
        if len(mask_file)>0:
            mask_file=mask_file[0]
            mask_file=os.path.join(image_file_path,mask_file)
        else:
            mask_file='mask'
        
        img_HE = load_image(HE_file)
        if os.path.exists(mask_file):
            mask=load_image(mask_file)
            if len(mask.shape)>2:
                mask=mask[:,:,0]
            mask=(mask>0).astype('float')
        else:
            if image_pre:
                mask=get_HE_mask_3d(HE_file,img_HE,image_file_path)
            elif mask_sam:
                import torch
                from segment_anything import sam_model_registry, SamAutomaticMaskGenerator, SamPredictor
                sam_checkpoint = "/histor/zhao/pangkun/JupyterLab/Task/github/segment-anything/checkpoint/sam_vit_h_4b8939.pth"
                model_type = "vit_h"
                device = (
                    "cuda"
                    if torch.cuda.is_available()
                    else "mps"
                    if torch.backends.mps.is_available()
                    else "cpu"
                )
                sam = sam_model_registry[model_type](checkpoint=sam_checkpoint)
                sam.to(device=device)
                mask_generator = SamAutomaticMaskGenerator(sam)
            else:
                mask=get_HE_mask(img_HE)
        
        if only_HE:
            line_point_r1c1=(line_point_r1c1[0],line_point_r1c1[1])
            line_point_r1c70=(line_point_r1c70[0],line_point_r1c70[1])
            line_point_r70c1=(line_point_r70c1[0],line_point_r70c1[1])
            
            #img_HE[mask==0]=255
            img = img_HE
            x, y = img.shape[0:2]

            tissue_hires_scalef = round(2000/x,8)
            tissue_hires_image = cv2.resize(img, (int(y * tissue_hires_scalef), 2000))
            if mask_sam:
                masks = mask_generator.generate(tissue_hires_image)
                sorted_anns = sorted(masks, key=(lambda x: x['area']), reverse=True)
                tissue_mask_image=sorted_anns[0]['segmentation'].astype(float)
            else:
                tissue_mask_image = cv2.resize(mask, (int(y * tissue_hires_scalef), 2000))
            tissue_hires_spot=np.zeros((2,2))

            tissue_lowres_scalef = round(600/x,8)
            tissue_lowres_image = cv2.resize(img, (int(y * tissue_lowres_scalef), 600))

            HE_row, HE_col = img.shape[0:2]
            del img_HE
            del img
            
        else:
            img_Spot_marker = load_image(spot_marker_file)
            img_spot = load_image(spot_file)
            #图像缩放
            d=dis_two_point(HE_point1,HE_point2)/dis_two_point(Spot_point1,Spot_point2)
            img_Spot_marker=cv2.resize(img_Spot_marker,[int(img_Spot_marker.shape[1]*d),int(img_Spot_marker.shape[0]*d)])
            img_spot=cv2.resize(img_spot,[int(img_spot.shape[1]*d),int(img_spot.shape[0]*d)])
            Spot_point1=(Spot_point1[0]*d,Spot_point1[1]*d)
            Spot_point2=(Spot_point2[0]*d,Spot_point2[1]*d)

            rows,cols,z = img_Spot_marker.shape
            # x轴平移200，y轴平移100, 2*3矩阵 
            r=-(Spot_point1[0]-HE_point1[0]+Spot_point2[0]-HE_point2[0])/2#-260#
            c=-(Spot_point1[1]-HE_point1[1]+Spot_point2[1]-HE_point2[1])/2#-221
            M = np.float32([[1, 0, c], [0, 1, r]]) #第1个 左右平移  第2个 上下平移
            # 用仿射变换实现平移
            img_Spot_marker = cv2.warpAffine(img_Spot_marker, M, (cols, rows)) #cols, row 正常的 row col
            img_spot = cv2.warpAffine(img_spot, M, (cols, rows)) #cols, row 正常的 row col

            #img_HE[mask==0]=255
            img = img_HE
            x, y = img.shape[0:2]

            tissue_hires_scalef = round(2000/x,8)
            tissue_hires_image = cv2.resize(img, (int(y * tissue_hires_scalef), 2000))
            if mask_sam:
                masks = mask_generator.generate(tissue_hires_image)
                sorted_anns = sorted(masks, key=(lambda x: x['area']), reverse=True)
                tissue_mask_image=sorted_anns[0]['segmentation'].astype(float)
            else:
                tissue_mask_image = cv2.resize(mask, (int(y * tissue_hires_scalef), 2000))

            xx, yy = img_spot.shape[0:2]
            tissue_hires_spot = cv2.resize(img_spot, (int(yy * tissue_hires_scalef), int(xx * tissue_hires_scalef)))

            tissue_lowres_scalef = round(600/x,8)
            tissue_lowres_image = cv2.resize(img, (int(y * tissue_lowres_scalef), 600))
            #tissue_mask_image = cv2.resize(mask, (int(y * tissue_lowres_scalef), 600))

            HE_row, HE_col = img.shape[0:2]
            del img_spot
            del img_Spot_marker
            del img_HE
            del img

            line_point_r1c1=(line_point_r1c1[0]*d+r,line_point_r1c1[1]*d+c)
            line_point_r1c70=(line_point_r1c70[0]*d+r,line_point_r1c70[1]*d+c)
            line_point_r70c1=(line_point_r70c1[0]*d+r,line_point_r70c1[1]*d+c)

    gene_info=pd.read_csv(geneinfo_path,sep='\t',index_col='Geneid')
    if chip_type in ['T3','T9']:
        X_raw=mmread(gene_file).astype(int)
        df_barcodes=pd.read_csv(barcodes_file,sep='\t',header=None)
        df_features=pd.read_csv(features_file,sep='\t',header=None)

        Spot_coords=get_chip_coord(barcode_num,reg,chip_type,Barcode_file_path)
        Barcode_path=Barcode_file_path+chip_type+'-spatial_barcode'+str(barcode_num)+'.txt'
        if img_file:
            pixel_coords,dis=get_piexl_coods(line_point_r1c1,line_point_r70c1,line_point_r1c70,barcode_num,reg,chip_type)
            pixel_coords=pd.DataFrame(pixel_coords,index=Spot_coords.index,columns=['col','row'])#如果需要导入HE
            
        Barcode=pd.read_csv(Barcode_path,sep='\t',header=None)
        Barcode.index=list(Barcode[1])
        Barcode=Barcode.loc[list(df_barcodes[0]),:]
        Spot_coords=Spot_coords.loc[list(Barcode[0]),:]
        gene_info=gene_info.loc[list(df_features[0]),:]
        if img_file:
            pixel_coords=pixel_coords.loc[list(Barcode[0]),:]
        
    elif chip_type in ['M9']:
        X_raw=mmread(gene_file).astype(int).toarray()
        df_barcodes=pd.read_csv(barcodes_file,sep='\t',header=None)
        df_features=pd.read_csv(features_file,sep='\t',header=None)
        sel_gene=np.array((X_raw>0).sum(1))>10
        df_features=df_features.loc[sel_gene,:]

        Barcode_path=Barcode_file_path+chip_type+'_ST_info_barcode_chip1_C18.csv'
        Barcode=pd.read_csv(Barcode_path,index_col=0)
        Barcode_uniq=Barcode.drop_duplicates(['coord'])
        Barcode=Barcode.sort_values('coord')

        df_barcodes[1]=list(range(df_barcodes.shape[0]))
        df_barcodes.index=list(df_barcodes[0])
        df_barcodes=df_barcodes.loc[list(Barcode.seq_ABC),:]
        
        print('Start merge-barcode-C')
        
        X_raw=X_raw.T[:,sel_gene]
        X_raw=X_raw[list(df_barcodes[1]),:]
        X_raw=sp.csr_matrix(X_raw[list(range(0,X_raw.shape[0],2)),:]+X_raw[list(range(1,X_raw.shape[0],2)),:]).T
        print('End merge-barcode-C')

        if Velocyto:
            x_spliced=x_spliced.A[:,sel_gene]
            x_spliced=x_spliced[list(df_barcodes[1]),:]
            x_spliced=sp.csr_matrix(x_spliced[list(range(0,x_spliced.shape[0],2)),:]+x_spliced[list(range(1,x_spliced.shape[0],2)),:])
            
            x_unspliced=x_unspliced.A[:,sel_gene]
            x_unspliced=x_unspliced[list(df_barcodes[1]),:]
            x_unspliced=sp.csr_matrix(x_unspliced[list(range(0,x_unspliced.shape[0],2)),:]+x_unspliced[list(range(1,x_unspliced.shape[0],2)),:])
        Barcode=Barcode.drop_duplicates(['coord'])
        if img_file:
            if add_M9_15 or add_M9_20:
                pixel_coords,pixel_coords_all,spot_coords_all,dis=get_piexl_coods_Mchip(line_point_r1c1,line_point_r70c1,line_point_r1c70,channels_num,barcode_num,chip_type,Barcode_file_path,add_M9_15,add_M9_20)
                pixel_coords=pd.DataFrame(pixel_coords,index=list(Barcode_uniq['coord']),columns=['col','row'])
            else:
                pixel_coords,dis=get_piexl_coods_Mchip(line_point_r1c1,line_point_r70c1,line_point_r1c70,channels_num,barcode_num,chip_type,Barcode_file_path,add_M9_15,add_M9_20)
                pixel_coords=pd.DataFrame(pixel_coords,index=list(Barcode_uniq['coord']),columns=['col','row'])
        Spot_coords=Barcode.loc[:,['col','row']]
        Spot_coords.index=list(Barcode['coord'])
        Barcode.index=list(Barcode['coord'])
        if img_file:
            pixel_coords=pixel_coords.loc[list(Barcode['coord']),:]

    obs = pd.DataFrame()
    obs.index=Barcode.index
    obs['sample'] = [sample for i in obs.index]
    if chip_type in ['T3','T9']:
        obs['reg'] = [reg for i in obs.index]

    obs['Spot_col'] = list(Spot_coords['col'])
    obs['Spot_row'] = list(Spot_coords['row'])
    if add_M9_15 or add_M9_20:
        df_tmp=pd.DataFrame(np.zeros((spot_coords_all.shape[0],obs.shape[1])))
        df_tmp.iloc[:,0:3]=np.hstack([np.array(['add' for i in range(spot_coords_all.shape[0])]).reshape((spot_coords_all.shape[0],1)),spot_coords_all])
        df_tmp.index=[ str(spot_coords_all[i,0])+'x'+str(spot_coords_all[i,1]) for i in range(spot_coords_all.shape[0])]
        df_tmp.columns=obs.columns
        df_tmp['Spot_col']=df_tmp['Spot_col'].astype('int')
        df_tmp['Spot_row']=df_tmp['Spot_row'].astype('int')
        obs=pd.concat( [obs,df_tmp])
    var=pd.DataFrame()
    var.index=list(gene_info.index)
    var['Symbol'] = list(gene_info['Symbol'])
    var['Gene_type'] = list(gene_info['type'])
    if  chip_type in ['M9']:
        var= var.loc[list(df_features[0]),:]
    
    Spot_spatial=np.array(Spot_coords[['col','row']])
    
    if img_file:
        Pixel_spatial=np.array(pixel_coords[['col','row']])
        if add_M9_15 or add_M9_20:
            Pixel_spatial=np.vstack([Pixel_spatial,pixel_coords_all])
            Spot_spatial=np.vstack([Spot_spatial,spot_coords_all])
            X_raw=sp.csr_matrix(np.vstack((X_raw.T.A, np.zeros((spot_coords_all.shape[0],var.shape[0])).astype(int))).T)
            if Velocyto:
                x_spliced=sp.csr_matrix(np.vstack((x_spliced.A, np.zeros((spot_coords_all.shape[0],var.shape[0])).astype(int))))
                x_unspliced=sp.csr_matrix(np.vstack((x_unspliced.A, np.zeros((spot_coords_all.shape[0],var.shape[0])).astype(int))))
        if mask_sam:
            mask=cv2.resize(tissue_mask_image, (HE_col,HE_row))
            tissue_in=[str(int(mask[i[1],i[0]])) for i in Pixel_spatial]
        else:
            tissue_in=[str(int(mask[i[1],i[0]])) for i in Pixel_spatial]
        
        obs['in_tissue']=tissue_in
        adata = AnnData(X_raw.T, obs=obs, var=var,obsm={"spatial": Pixel_spatial},
                       uns={'spatial':{sample:{'images':{'hires':tissue_hires_image,
                                                         'lowres':tissue_lowres_image,
                                                         'hires_mask':tissue_mask_image,
                                                          'hires_spot':tissue_hires_spot
                                                             },
                                                'scalefactors':{'tissue_hires_scalef': tissue_hires_scalef,
                                                               'tissue_lowres_scalef': tissue_lowres_scalef,
                                                               'tissue_hires_spot_scalef': tissue_hires_scalef,
                                                                #'fiducial_diameter_fullres': dis,
                                                                'spot_diameter_fullres': round(dis,8),
                                                                'microns_pixels': round(res_um/dis,8) }
                                               }}})
        if Velocyto:
            adata.layers={'spliced':x_spliced, 'unspliced':x_unspliced}
        adata.uns['HE_row']=HE_row
        adata.uns['HE_col']=HE_col
        if chip_type in ['T3','T9']:
            adata.X=sp.csr_matrix(adata.X.toarray())

        adata.var_names_make_unique()

        if species == 'human':
            adata.var["mt"] = adata.var.Symbol.str.startswith("MT-")
            adata.var['rp'] = adata.var.Gene_type.str.startswith("rRNA")
            adata.var["hb"] = [True if i.startswith('HB') and not i.startswith('HBP') else False for i in adata.var.Symbol]
        elif species == 'mouse':
            adata.var["mt"] = adata.var.Symbol.str.startswith("mt-")
            adata.var['rp'] = adata.var.Gene_type.str.startswith("rRNA")
            adata.var["hb"] = [True if i.startswith('hb') and not i.startswith('hbp') else False for i in adata.var.Symbol]
        
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","rp","hb"], inplace=True)
        print(f"#Spot in tissue: {sum(adata.obs['in_tissue']=='1')}")
        print(f"#Spot in tissue ratio: {int(100*sum(adata.obs['in_tissue']=='1')/(channels_num*channels_num))}%")
        print('Median UMIs:',round(np.median(adata[adata.obs["in_tissue"] == '1' ].obs['total_counts']),0))
        print('Median Genes:',round(np.median(adata[adata.obs["in_tissue"] == '1' ].obs['n_genes_by_counts']),0))
        print('In tissue(%):',round(100-100*np.sum(adata.obs['total_counts'][adata.obs['in_tissue']=='0'])/np.sum(adata.obs['total_counts']),0))
    else:
        adata = AnnData(X_raw.T, obs=obs, var=var,obsm={"spatial": Spot_spatial},
        uns={'spatial':{sample:{'image':{'hires':np.zeros([2,2])}}}})
        if Velocyto:
            adata.layers={'spliced':x_spliced, 'unspliced':x_unspliced}
        
        if chip_type in ['T3','T9']:
            adata.X=sp.csr_matrix(adata.X.toarray())
        adata.var_names_make_unique()
        #adata.layers['raw']=copy.deepcopy(adata.X.astype(int))
        if species == 'human':
            adata.var["mt"] = adata.var.Symbol.str.startswith("MT-")
            adata.var['rp'] = adata.var.Gene_type.str.startswith("rRNA")
            adata.var["hb"] = [True if i.startswith('HB') and not i.startswith('HBP') else False for i in adata.var.Symbol]
        elif species == 'mouse':
            adata.var["mt"] = adata.var.Symbol.str.startswith("mt-")
            adata.var['rp'] = adata.var.Gene_type.str.startswith("rRNA")
            adata.var["hb"] = [True if i.startswith('hb') and not i.startswith('hbp') else False for i in adata.var.Symbol]
        
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","rp","hb"], inplace=True)
        print('No Image')
    end_time = time.time()
    print("Reading runtime：", round(end_time - start_time,0), "s")
    print('-')
    return adata

def add_velo(adata_raw):
    import scvelo as scv
    adata_raw=adata_raw.copy()
    loom_file='st_data/'+list(adata_raw.uns['spatial'].keys())[0]+'.loom'
    ldata = scv.read(loom_file, cache=False)
    ldata.obs.index=[ i.split(':')[1][0:16]  for i in ldata.obs.index]
    ldata.var.index=list(ldata.var['Accession'])
    adata_velo = scv.utils.merge(adata_raw, ldata)
    return adata_velo



def read_10X_visium(file_path,count_file):
    import scanpy as sc
    import copy
    adata=sc.read_visium(file_path,count_file=count_file)
    adata.obsm['spatial']=adata.obsm['spatial'].astype(int)
    adata.obs['in_tissue']=adata.obs['in_tissue'].astype(str)
    sample=list(adata.uns['spatial'].keys())[0]
    adata.var_names_make_unique()
    adata.layers['raw']=copy.deepcopy(adata.X)
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["hb"] = adata.var_names.str.startswith("Hb-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","hb"], inplace=True)
    print(f'{sample}:')
    print('Median total_counts',np.median(adata[adata.obs["in_tissue"] == '1'  ].obs['total_counts']))
    print('Median n_genes_by_counts',np.median(adata[adata.obs["in_tissue"] == '1'  ].obs['n_genes_by_counts']))
    print('-')
    return adata

def get_adata_STARsolo_10x(sample,EM=True,Velocyto=False,soloFeatures='Gene',raw=True,species='human',s_id=False,s=None,
                        file_path=None,image_file_path=None,Barcode_file_path=None,geneinfo_path=None):
    import time
    start_time = time.time()
    import cv2 
    import scanpy as sc
    import os
    import numpy as np
    import pandas as pd
    from scipy.io import mmread
    from anndata import AnnData
    import scipy.sparse as sp
    print(f'sample: {sample}')
    if s_id:
        sample_name=sample[:sample.rfind(s)]
    scale_file=Barcode_file_path+sample_name+'-scalefactors_json.json'
    hires_file=Barcode_file_path+sample_name+'-tissue_hires_image.png'
    lowres_file=Barcode_file_path+sample_name+'-tissue_lowres_image.png'
    import json
    with open(scale_file, 'r', encoding='utf-8') as file:  
        scale_data = json.load(file) 
    tissue_hires_image=load_image(hires_file)
    tissue_lowres_image=load_image(lowres_file)
    if Velocyto:
        print('Reading Velocyto')         
        if raw:
            #x_ambiguous=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/raw/ambiguous.mtx.gz').astype(int).toarray().T)
            x_spliced=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/raw/spliced.mtx.gz').astype(int).toarray().T)
            x_unspliced=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/raw/unspliced.mtx.gz').astype(int).toarray().T)
        else:
            #x_ambiguous=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/filtered/ambiguous.mtx.gz').astype(int).toarray().T)
            x_spliced=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/filtered/spliced.mtx.gz').astype(int).toarray().T)
            x_unspliced=sp.csr_matrix(mmread(file_path+sample+'_Solo.out/Velocyto/filtered/unspliced.mtx.gz').astype(int).toarray().T)
    else:
        print('No Velocyto')         

    if soloFeatures=='Gene' or soloFeatures=='GeneFull':
        if raw:
            if EM:
                gene_file=file_path+sample+'_Solo.out/'+soloFeatures+'/raw/UniqueAndMult-EM.mtx.gz'
            else:
                gene_file=file_path+sample+'_Solo.out/'+soloFeatures+'/raw/matrix.mtx.gz'
            barcodes_file=file_path+sample+'_Solo.out/'+soloFeatures+'/raw/barcodes.tsv'
            features_file=file_path+sample+'_Solo.out/'+soloFeatures+'/raw/features.tsv'
        else:
            gene_file=file_path+sample+'_Solo.out/'+soloFeatures+'/filtered/matrix.mtx.gz'
            barcodes_file=file_path+sample+'_Solo.out/'+soloFeatures+'/filtered/barcodes.tsv'
            features_file=file_path+sample+'_Solo.out/'+soloFeatures+'/filtered/features.tsv'
    else:
        print('No soloFeatures')
    geneinfo_path=geneinfo_path+'geneInfo.txt'
    gene_info=pd.read_csv(geneinfo_path,sep='\t',index_col='Geneid')
    X_raw=mmread(gene_file).astype(int)
    df_barcodes=pd.read_csv(barcodes_file,sep='\t',header=None)
    df_features=pd.read_csv(features_file,sep='\t',header=None)
    Barcode_path=Barcode_file_path+sample_name+'-tissue_positions_list.csv'
    Barcode=pd.read_csv(Barcode_path)
    Barcode.index=[i[:-2] for i in Barcode['barcode']]
    obs=Barcode.loc[list(df_barcodes[0]),:]
    obs['in_tissue'] = [str(i) for i in obs['in_tissue']]
    pixel_spatial=np.array(obs[['pxl_col_in_fullres','pxl_row_in_fullres']])
    gene_info=gene_info.loc[list(df_features[0]),:]

    var=pd.DataFrame()
    var.index=list(gene_info.index)
    var['Symbol'] = list(gene_info['Symbol'])
    var['Gene_type'] = list(gene_info['type'])

    adata = AnnData(X_raw.T, obs=obs, var=var,obsm={"spatial": pixel_spatial},
                   uns={'spatial':{sample:{'images':{'hires':tissue_hires_image,
                                                     'lowres':tissue_lowres_image
                                                         },
                                            'scalefactors':scale_data
                                           }}})
    if Velocyto:
        adata.layers={'spliced':x_spliced, 'unspliced':x_unspliced}
    adata.X=sp.csr_matrix(adata.X.toarray())
    adata.var_names_make_unique()
    if species == 'human':
        adata.var["mt"] = adata.var.Symbol.str.startswith("MT-")
        adata.var['rp'] = adata.var.Gene_type.str.startswith("rRNA")
        adata.var["hb"] = [True if i.startswith('HB') and not i.startswith('HBP') else False for i in adata.var.Symbol]
    elif species == 'mouse':
        adata.var["mt"] = adata.var.Symbol.str.startswith("mt-")
        adata.var['rp'] = adata.var.Gene_type.str.startswith("rRNA")
        adata.var["hb"] = [True if i.startswith('hb') and not i.startswith('hbp') else False for i in adata.var.Symbol]
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","rp","hb"], inplace=True)
    print(f"#Spot in tissue: {sum(adata.obs['in_tissue']=='1')}")
    print(f"#Spot in tissue ratio: {int(100*sum(adata.obs['in_tissue']=='1')/4992)}%")
    print('Median UMIs:',round(np.median(adata[adata.obs["in_tissue"] == '1' ].obs['total_counts']),0))
    print('Median Genes:',round(np.median(adata[adata.obs["in_tissue"] == '1' ].obs['n_genes_by_counts']),0))
    print('In tissue(%):',round(100-100*np.sum(adata.obs['total_counts'][adata.obs['in_tissue']=='0'])/np.sum(adata.obs['total_counts']),0))
    end_time = time.time()
    print("Reading runtime：", round(end_time - start_time,0), "s")
    print('-')
    return adata
def get_adata_st_100um(sample,file_path,res_um,species):
    import time
    import os
    start_time = time.time()
    import cv2 
    import scanpy as sc
    import numpy as np
    import pandas as pd
    from anndata import AnnData
    import scipy.sparse as sp
    df_count=pd.read_csv(file_path+'/count-matrices/'+sample+'.tsv.gz',sep='\t',index_col=0)
    #df_label=pd.read_csv(file_path+'/meta/'+sample+'_labeled_coordinates.tsv',sep='\t',index_col=0)
    df_selection=pd.read_csv(file_path+'/spot-selections/'+sample+'_selection.tsv.gz',sep='\t')
    df_selection.index=[ str(df_selection['x'][i])+'x'+ str(df_selection['y'][i]) for i in range(df_selection.shape[0])]
    df_selection=df_selection.loc[df_count.index,:]
    img = cv2.imread(file_path+'/images/HE/'+sample+'.jpg')[:,:,::-1]
    x, y = img.shape[0:2]
    tissue_hires_scalef = round(2000/x,8)
    tissue_hires_image = cv2.resize(img, (int(y * tissue_hires_scalef), 2000))
    if os.path.exists(file_path+'/images/annotation/'+sample+'.jpg'):
        img_HE_anno = cv2.imread(file_path+'/images/annotation/'+sample+'.jpg')[:,:,::-1]
        #anno_x=int(2000/(img_HE_anno.shape[1]/img.shape[1]))
        ss=img_HE_anno.shape[1]/img.shape[1]
        x, y = img_HE_anno.shape[0:2]
        tissue_lowres_scalef = round(2000/x,8)
        tissue_lowres_image = cv2.resize(img_HE_anno, (int(y * tissue_lowres_scalef), 2000))
        tissue_lowres_scalef=tissue_lowres_scalef*ss
        del img_HE_anno
    else:
        tissue_lowres_image=tissue_hires_image
        tissue_lowres_scalef=tissue_hires_scalef
        
    HE_row, HE_col = img.shape[0:2]
    
    del img
    
    dis = (max(df_selection['pixel_y'])-min(df_selection['pixel_y']))/(max(df_selection['y'])-min(df_selection['y']))/2
    var=pd.DataFrame()
    var.index=list(df_count.columns)
    Pixel_spatial=np.array(df_selection[['pixel_x','pixel_y']])
    
    df_selection["in_tissue"] = '1'
    
    adata = AnnData(df_count, obs=df_selection, var=var,obsm={"spatial": Pixel_spatial},
                   uns={'spatial':{sample:{'images':{'hires':tissue_hires_image,
                                                     'lowres':tissue_lowres_image,
                                                         },
                                            'scalefactors':{'tissue_hires_scalef': tissue_hires_scalef,
                                                           'tissue_lowres_scalef': tissue_lowres_scalef,
                                                            'spot_diameter_fullres': round(dis,8),
                                                            'microns_pixels': round(res_um/dis,8) }
                                           }}})
    adata.X=sp.csr_matrix(adata.X)
    if species == 'human':
        adata.var["mt"] = adata.var.index.str.startswith("MT-")
        adata.var["hb"] = [True if i.startswith('HB') and not i.startswith('HBP') else False for i in adata.var.index]
    elif species == 'mouse':
        adata.var["mt"] = adata.var.index.str.startswith("mt-")
        adata.var["hb"] = [True if i.startswith('hb') and not i.startswith('hbp') else False for i in adata.var.index]
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","hb"], inplace=True)
    
    adata.uns['HE_row']=HE_row
    adata.uns['HE_col']=HE_col
    
    print(f'{sample}')
    print('Median UMIs:',round(np.median(adata.obs['total_counts']),0))
    print('Median Genes:',round(np.median(adata.obs['n_genes_by_counts']),0))
    end_time = time.time()
    print("Reading runtime：", round(end_time - start_time,0), "s")
    print('-')
    return adata
    
def adata_filter_norm_process(adata_raw,min_cells_num=10,min_genes_num=100,mt_num=20,n_top_genes=12000,use_variable=True,
                              in_tissue=True,image=True,n_clusters=8,scale=False): 
    import scanpy as sc
    adata=adata_raw.copy()
    if in_tissue:
        adata = adata[adata.obs["in_tissue"] == '1' ]
    if image:
        adata.uns['crop_coord']= [min(adata.obsm['spatial'][:,0])-100, max(adata.obsm['spatial'][:,0])+100,
               min(adata.obsm['spatial'][:,1])-100, max(adata.obsm['spatial'][:,1])+100]

    else:
         adata.uns['crop_coord']= [min(adata.obsm['spatial'][:,0])-3, max(adata.obsm['spatial'][:,0])+3,
               min(adata.obsm['spatial'][:,1])-3, max(adata.obsm['spatial'][:,1])+3]
    if scale:
        sc.pp.filter_cells(adata, min_genes=min_genes_num)
        sc.pp.filter_genes(adata, min_cells=min_cells_num)
        adata = adata[adata.obs["pct_counts_mt"] < mt_num]
        
        adata.layers['raw'] = adata.X.copy()
        sc.pp.normalize_total(adata,target_sum=1e4)
        adata.layers['norm'] = adata.X.copy()
        sc.pp.log1p(adata)
        adata.raw = adata.copy()
        sc.pp.highly_variable_genes(adata,n_top_genes=n_top_genes) 
        adata=adata[:,adata.var.highly_variable]#adata.raw.to_adata()
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)
        
    else:
        sc.pp.filter_cells(adata, min_genes=min_genes_num)
        sc.pp.filter_genes(adata, min_cells=min_cells_num)
        adata = adata[adata.obs["pct_counts_mt"] < mt_num]
        adata.layers['raw'] = adata.X.copy()
        sc.pp.normalize_total(adata,target_sum=1e4)
        adata.raw = adata.copy()
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata,n_top_genes=n_top_genes) 

    sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=use_variable)
    sc.pp.neighbors(adata,  use_rep='X_pca', n_neighbors=50)
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    sc.tl.leiden(adata)
    sc.tl.louvain(adata, resolution=1)
    X_pca = adata.obsm['X_pca'] 
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X_pca) 
    adata.obs['kmeans'] = kmeans.labels_.astype(str)
    
    print('Sample:  '+list(adata.uns['spatial'].keys())[0])
    print(f"#genes after  filter: {adata.n_vars}")
    print(f"#Spot after  filter: {adata.n_obs}")
    print(f"#Median UMIs : {round( np.median(adata.obs['total_counts']),0)}")
    print(f"#Median Genes : {round(np.median(adata.obs['n_genes_by_counts']),0)}")
    print(f"#Mean UMIs : {round( np.mean(adata.obs['total_counts']),0)}")
    print(f"#Mean Genes : {round(np.mean(adata.obs['n_genes_by_counts']),0)}")
    return adata
