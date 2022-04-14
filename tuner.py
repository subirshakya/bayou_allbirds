import subprocess
import os
import glob
import csv

os.chdir("C:/Subir/Projects/Avonet/")

def run_r(group, alpha, beta, sigma, theta):
    files = glob.glob('C:/Subir/Projects/Avonet/runs/test/*')
    #print(files)
    for file in files:
        os.remove(file)
    myoutput = open('test.txt', 'w')
    subprocess.call (["Rscript", "ou_test.R", str(group), str(alpha), str(sigma), str(beta), str(theta)], stdout = myoutput, stderr=subprocess.DEVNULL)
    myoutput.close()

def tuner(alpha, a, beta, b, sigma, s, theta, t):
    if a < 0.2:
        alpha -= 1
    elif a > 0.4:
        alpha += 1
    if alpha <= 0:
        alpha = 1

    if b < 0.2:
        beta -= 0.1
    elif b > 0.4:
        beta += 0.1
    if beta <= 0:
        beta = 0.1      
        
    if s < 0.2:
        sigma -= 0.2
    elif s > 0.4:
        sigma += 0.2
    if sigma <= 0:
        sigma = 0.5        

    if t < 0.2:
        theta -= 0.5
    elif t > 0.4:
        theta += 0.5
    if theta <= 0:
        theta = 0.5        
    return(alpha, beta, sigma, theta)

def get_vals():
    reader = open('test.txt', 'r')
    for row in reader.readlines():
        last_row = row.split()
    reader.close()

    return(float(last_row[8]), float(last_row[9]), float(last_row[10]), float(last_row[11])) 

group = 18
alpha = 6
beta = 0.3
sigma = 1
theta = 3
run_r(group, alpha, beta, sigma, theta)
a,b,s,t = get_vals()

while (a < 0.2 or a > 0.4) or (b < 0.2 or b > 0.4) or (s < 0.2 or s > 0.4) or (t < 0.2 or t > 0.4):
    print(a,b,s,t,alpha,beta,sigma,theta)
    alpha, beta, sigma, theta = tuner(alpha, a, beta, b, sigma, s, theta, t)
    run_r(group, alpha, beta, sigma, theta)
    a,b,s,t = get_vals()

print("Final = ", alpha,beta,sigma,theta)
