# download file from the site 
from time import timezone
from datetime import datetime, timezone


def get_file(file_num):
    import requests
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.txt'
    response = requests.get(url)
    with open(f"./coords{file_num}.txt", "wb") as file :
        file.write(response.content)
 
# get state vector from the file 
def get_state_vector(y,m,d,h,mins,file_num):

    # if mins % 4 != 0 :
    #     mins -= (mins%4)
    if h <= 9 :
        h = f"0{h}"
    if m <= 9 :
        m = f"0{m}"
    if mins <= 9 :
        mins = f"0{mins}"

    date = f'{y}-{m}-{d}T{h}:{mins}'
    
    # print(date)
    with open(rf'./coords{file_num}.txt', 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(date) != -1:
                idx  = lines.index(line)
                line = line.replace("\n", "")
                line_arr = line.split(" ")
                r = [float(i) for i in line_arr[1:4]]
                v = [float(i) for i in line_arr[4:]]
                
                return idx , r, v 

    return 0

def get_by_idx(file_num , idx):

    with open(rf'./coords{file_num}.txt', 'r') as fp:
        lines = fp.readlines()
        line = lines[idx]
        line = line.replace("\n", "")
        line_arr = line.split(" ")
        r = [float(i) for i in line_arr[1:4]]
        v = [float(i) for i in line_arr[4:]]        
        return  r, v 

def get_state_vector_now(file_num):
    now = datetime.now(timezone.utc) #current date and time
    year = int(now.strftime("%Y"))
    month = int(now.strftime("%m"))
    day = int(now.strftime("%d"))
    time = now.strftime("%H:%M").split(":")
    hour = int(time[0])
    mins = int(time[1])
    result = get_state_vector(year,month,day,hour,mins,file_num)

    if result :
        return result
    else:
        for i in range(6):
            mins = mins + 1
            if  get_state_vector(year,month,day,hour,mins,file_num) :
                return  get_state_vector(year,month,day,hour,mins,file_num)

