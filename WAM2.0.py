import os
import re
import numpy
import math
import random
import matplotlib.pyplot as plt
"""
上机实践1：真核基因剪接位点预测的WAM模型的编程实践
选择基因末端前3个位点与后6个位点作为donor位点
选择基因前端前3个位点与后6个位点作为acceptor位点
在序列的第一位到第lengh-9位上挑选长度为9的序列作为negative  sequence
首先根据training set 构建 P+/P-矩阵 8*16 AA AC AG AT.....
利用构建出来的P+/P-矩阵，根据s(x)公式对training set 中的donor位点进行计算，求出求平均值，作为判别donor/accep位点的阈值
对testing set中的数据进行甄别，通过s(x)值来统计TP、FP、FN来计算sn，sp。通过改变s(x)的范围来得到不同的sn，sp的值，从而做出sn-sp曲线
"""
###构建赋值函数
def get_value(S):
    """
    用于定位到所要构建的P+,P-矩阵的特定的行
    利用 4*get_value(S)+get_value(s+1)这个函数来定位AA，AC,AG、、、所在的列
    :param S:某个字符
    :return:
    """
    if S == "A":
        a = 0
    if S == "C":
        a = 1
    if S == "G":
        a = 2
    if S == "T":
        a = 3
    return a

###构建p矩阵函数
def p_matrix(sequence):
    """
    根据字典sequence构造p_matrix,并返回单核苷酸矩阵与双核苷酸矩阵
    :param sequence:
    :return:
    """
    snp_count = numpy.zeros((9,4))
    dnp_count = numpy.zeros((8,16))
    dnp_p_matrix = numpy.zeros((8,16))
    for i in range(len(sequence)):
        for n in range(len(sequence[i])):
            if sequence[i][n] == "A":
                snp_count[n][0] += 1
            if sequence[i][n] == "C":
                snp_count[n][1] += 1
            if sequence[i][n] == "G":
                snp_count[n][2] += 1
            if sequence[i][n] == "T":
                snp_count[n][3] += 1
            if n < 8:
                b = 4 * get_value(sequence[i][n]) + get_value(sequence[i][n + 1])
                dnp_count[n][b] += 1
    for i in range(len(dnp_count)):
        for j in range(len(dnp_count[i])):
            if snp_count[i][j // 4] / max(snp_count[i][0], snp_count[i][1], snp_count[i][2], snp_count[i][3]) < 0.1 and i != 2:
                dnp_p_matrix[i][j] = 0

            else:
                dnp_p_matrix[i][j] = (dnp_count[i][j] / snp_count[i][j // 4]).round(4)
    return snp_count, dnp_p_matrix

def calculation_of_sx(Ppos, Ppos_1, Pneg, Pneg_1, s):
    """
    计算每个片段的评分S(x)
    :param Ppos:p+矩阵，双核苷酸矩阵
    :param Ppos_1:第一位单核苷酸的概率值
    :param Pneg:p-矩阵，双核苷酸矩阵
    :param Pneg_1:第一位核苷酸的概率值
    :param s: 存储有多个片段的字典
    :return:
    """
    Sx = 0
    for n in range(len(s)):
        if s[n] != "A" and s[n] != "C" and s[n] != "G" and s[n] != "T":
            Sx = 1000
            return Sx
    for i in range(len(s)):
        b = get_value(s[i - 1]) * 4 + get_value(s[i])
        if i == 0:
            Sx = math.log(Ppos_1[get_value(s[0])] / Pneg_1[get_value(s[0])]) + Sx
            #计算第一位 log(1,P+)/log(1,P-)
        elif Pneg[i - 1][b] == 0 or Ppos[i - 1][b] == 0:
            #若矩阵中该位置概率为0则+math.log(1)
            Sx = Sx + math.log(1)
        else:
            Sx = math.log(Ppos[i - 1][b] / Pneg[i - 1][b]) + Sx
    return Sx

def jugde(s):
    """
    判断此片段是否含有N碱基，若全为ACGT则返回0，否则返回1
    """
    a = 0
    for i in range(len(s)):
        if s[i] != "A" and s[i] != "C" and s[i] != "G" and s[i] != "T":
            a = 1
    return a

train_donor_num = 0
train_random_num = 0
train_acceptor_num = 0
Donor = {}
Acceptor = {}
train_random = {}
TP = 0
FP = 0
FN = 0


root_dir = r"D:\JetBrains\shujuwajue\Training Set"
if True:
    for file in os.listdir(root_dir):
        file_name = root_dir + "\\" + file
        filein = open(file_name, "r")
        doner_pos = []
        acceptor_pos = []
        doner_seq = []
        acceptor_seq = []
        random_seq = []
        xulie = str()

        for line in filein:
            if line[0:3] == "CDS":
                pos = re.sub("\D", " ", line)
                # print(line)
                pos = pos.split()
                # print(pos)
                for i in range(len(pos)):
                    if i % 2 != 0:
                        doner_pos.append(pos[i])
                        # 得到doner 位点的端点序号 并保存在一个元组中
                    if i % 2 == 0:
                        acceptor_pos.append(pos[i])
                        # 得到 acceptor 位点的端点序号 并保存在一个元组中
            if line[0] == "a" or line[0] == "t" or line[0] == "c" or line[0] == "g":
                xulie = xulie + line[0:-1]
            # 得到这段read的总序列

        for i in range(1,len(doner_pos)-1):
            for m in range(5):
                random_seq_num = random.randint(1, len(xulie) - 9)
                random_seq = xulie[random_seq_num:random_seq_num + 9].upper()
                jugde_random = jugde(random_seq)
                if jugde_random == 0:
                    train_random[train_random_num] = random_seq
                    train_random_num = train_random_num + 1
            """
            得到 一系列随机的 9nt长的片段 
            将它们存入train_random这个字典中
            其数量为train_random_num
            """
            doner_seq = xulie[eval(doner_pos[i]) - 3:eval(doner_pos[i]) + 6].upper()
            if doner_seq[1] != "N" and len(doner_seq) == 9:
                Donor[train_donor_num] = doner_seq
                train_donor_num = train_donor_num + 1
            """
            得到 training set 中所有的 doner 位点，截取其中前3个nt，后6个nt的片段
            将它们存入Donor这个字典中
            其数目为train_donor_num
            """
            acceptor_seq = xulie[eval(acceptor_pos[i]) - 3:eval(acceptor_pos[i]) + 6].upper()
            if acceptor_seq[1] != "N" and len(acceptor_seq) == 9:
                Acceptor[train_acceptor_num] = acceptor_seq
                train_acceptor_num = train_acceptor_num + 1

            """
            得到 training set 中所有的 acceptor位点，截取其中前3个nt，后6个nt的片段
            将它们存入Acceptor这个字典中
            其数目为train_acceptor_num
            """

# 得到Donor 的P+矩阵和第一行的概率值
Donor_count_zheng,Donor_Pzheng = p_matrix(Donor)
Donor_Pzheng_1 = Donor_count_zheng[0] / train_donor_num
# 得到Acceptor 的P+矩阵和第一行的概率值
Acceptor_count_zheng,Acceptor_Pzheng = p_matrix(Acceptor)
Acceptor_Pzheng_1 = Acceptor_count_zheng[0] / train_acceptor_num
# 得到 P-矩阵和第一行的概率值
count_fu,Pfu = p_matrix(train_random)
Pfu_1 = count_fu[0] / train_random_num

#计算出Donor位点的Sx的平均期望
Donor_Sx = 0
for i in range(len(Donor)):
    Donor_Sx = calculation_of_sx(Donor_Pzheng, Donor_Pzheng_1, Pfu, Pfu_1, Donor[i]) + Donor_Sx
average_Donor_Sx = Donor_Sx /train_donor_num
print("Donor位点的S(x)阈值: {}".format(average_Donor_Sx))

#计算出Acceptor位点的Sx的平均期望
Acceptor_Sx = 0
for i in range(len(Acceptor)):
    Acceptor_Sx = calculation_of_sx(Acceptor_Pzheng,Acceptor_Pzheng_1,Pfu,Pfu_1,Acceptor[i]) + Acceptor_Sx
average_Acceptor_Sx = Acceptor_Sx / train_acceptor_num
print("Acceptor位点的S(x)阈值: {}".format(average_Acceptor_Sx))

print("##############")
test_donor_num = 0
test_random_num = 0
test_acceptor_num = 0
Donor_test = {}
random_test = {}
Acceptor_test = {}
root_dir = r"D:\JetBrains\shujuwajue\Testing Set"

if True:
    for file in os.listdir(root_dir):
        file_name = root_dir + "\\" + file
        filein = open(file_name, "r")
        xulie_test = str()
        donor_test_pos = []
        acceptor_test_pos = []
        for line in filein:
            if line[0] == ">" and len(line) > 25:
                weidian = str(re.findall(r'[(](.*?)[)]', line))
                test_pos = re.sub("\D", " ", weidian)
                test_pos = test_pos.split()
                for i in range(1,len(test_pos)-1):
                    if i % 2 != 0:
                        donor_test_pos.append(test_pos[i])
                        #提取testing set 中的 donor位点的序号，并存储在donor_test_pos 这个元组中
                    if i % 2 == 0:
                        acceptor_test_pos.append(test_pos[i])
                        # 提取testing set 中的 acceptor位点的序号，并存储在acceptor_test_pos 这个元组中
            if line[0] == "A" or line[0] == "T" or line[0] == "C" or line[0] == "G":
                xulie_test = xulie_test + line[0:-1]
                #将序列存储到 xulie_test中
        for i in range(len(donor_test_pos)):
            #产生3倍数目的random test 片段，并将其存入 random_test中，其数目为test_random_num
            for n in range(5):
                random_seq_num = random.randint(1, len(xulie_test) - 9)
                random_seq = xulie_test[random_seq_num:random_seq_num + 9].upper()
                jugde_random = jugde(random_seq)
                if jugde_random == 0:
                    random_test[test_random_num] = random_seq
                    test_random_num = test_random_num + 1

            #选取Donor位点前3个nt，后6个nt的片段，并将其存入Donor_test这个字典中，其数目为 test_donor_num
            Donor_test_seq = xulie_test[eval(donor_test_pos[i]) - 3:eval(donor_test_pos[i]) + 6].upper()
            judge_test_Donor = jugde(Donor_test_seq)
            if judge_test_Donor == 0 and eval(donor_test_pos[i]) < (len(xulie_test) - 6):
                Donor_test[test_donor_num] = Donor_test_seq
                test_donor_num = test_donor_num + 1

            #选取Acceptor位点前3个nt，后6个nt的片段，并将其存入Acceptor_test这个字典中，其数目为 test_acceptor_num
            Acceptor_test_seq = xulie_test[eval(acceptor_test_pos[i]) - 3:eval(acceptor_test_pos[i]) + 6].upper()
            judge_test_Acceptor = jugde(Acceptor_test_seq)
            if judge_test_Acceptor == 0 and eval(acceptor_test_pos[i]) < (len(xulie_test) - 6):
                Acceptor_test[test_acceptor_num] = Acceptor_test_seq
                test_acceptor_num = test_acceptor_num +1


####Sn=TP/(TP+FN) TP表示预测的真实的Donor剪接位点,FN表示未能预测到的但仍是真实的Donor剪接位点，
####Sp=TP/(TP+FP) FP表示预测的假的Donor位点
C = 0
print("测试集中，donor位点的数目为: {}, acceptor位点的数目为: {}, 随机片段的数目为: {} ".format(test_donor_num, test_acceptor_num,test_random_num)) #jjjj是测试集中donor位点的数目，jjj是测试集中随机挑选的阴性数据的数目
D_XSn = numpy.zeros(100)
D_XSp = numpy.zeros(100)
aa = 100
for C in range(11, 111):
    TP = 0
    FP = 0
    FN = 0
    for i in range(test_random_num):
        Cfen = calculation_of_sx(Donor_Pzheng, Donor_Pzheng_1, Pfu, Pfu_1, random_test[i])
        if abs(Cfen - average_Donor_Sx) < C / 18:
        # if math.log(abs(Cfen-average_Donor_Sx)/average_Donor_Sx) < C / 50:
            FP = FP + 1 #统计假阳性数据个数

    for i in range(test_donor_num):
        # print(Donor_test[i])
        Cfen = calculation_of_sx(Donor_Pzheng, Donor_Pzheng_1, Pfu, Pfu_1, Donor_test[i])
        if abs(Cfen - average_Donor_Sx) < C / 18:
        # if math.log(abs(Cfen-average_Donor_Sx)/average_Donor_Sx) < C / 50:
            TP = TP + 1
    FN = test_donor_num - TP
    Sn = TP / (TP + FN)
    Sp = TP / (TP + FP)
    aa = aa - 1
    D_XSn[aa] = Sn
    D_XSp[aa] = Sp
print(D_XSn, D_XSp)
plt.figure("donor_site")
plt.plot(D_XSn[0:-1], D_XSp[0:-1], "-")
plt.xlabel("Sn")
plt.ylabel("Sp")
plt.show()

A_XSn = numpy.zeros(100)
A_XSp = numpy.zeros(100)
aa = 100
for C in range(11, 111):
    TP = 0
    FP = 0
    FN = 0
    for i in range(test_random_num):
        Cfen = calculation_of_sx(Acceptor_Pzheng, Acceptor_Pzheng_1, Pfu, Pfu_1, random_test[i])
        if abs(Cfen - average_Acceptor_Sx) < C / 18:
        # if math.log(abs(Cfen-average_Donor_Sx)/average_Donor_Sx) < C / 50:
            FP = FP + 1 #统计假阳性数据个数

    for i in range(test_acceptor_num):
        # print(Donor_test[i])
        Cfen = calculation_of_sx(Acceptor_Pzheng, Acceptor_Pzheng_1, Pfu, Pfu_1, Acceptor_test[i])
        if abs(Cfen - average_Acceptor_Sx) < C / 18:
        # if math.log(abs(Cfen-average_Donor_Sx)/average_Donor_Sx) < C / 50:
            TP = TP + 1
    FN = test_acceptor_num - TP
    Sn = TP / (TP + FN)
    Sp = TP / (TP + FP)
    aa = aa - 1
    A_XSn[aa] = Sn
    A_XSp[aa] = Sp
print(A_XSn, A_XSp)
plt.figure("acceptor_site")
plt.plot(A_XSn[0:-1], A_XSp[0:-1], "-")
plt.xlabel("Sn")
plt.ylabel("Sp")
plt.show()