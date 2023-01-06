name = ['2010BH115.py','2010JX16.py','2010NS36.py']
hv = ['16.5','16.4','15.4']
for i in range(3):
    namei = name[i]
    hvi = hv[i]
    f = open(namei, "r", encoding="utf-8")
    str1 = f.read()
    str2 = str1.replace('99.99',hvi)
    str3 = str2.replace("D_gs = 1","D_gs = 500 #")
    str4 = str3.replace("D_gss = [0,","D_gss = [1,8000]#")
    f.close()
    ff = open(namei, "w")
    ff.write(str4)
    ff.flush()
    ff.close()
