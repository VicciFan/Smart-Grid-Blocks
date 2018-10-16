import socket
import select

HOST = '2002:80b2:9c51::ffff:901'
PORT = 5004

s = socket.socket(socket.AF_INET6, socket.SOCK_DGRAM)
repeat = 0
total = 0
mean = 0

s.setblocking(0)

while repeat <= 59:
    i = 0
    while True:
        i = i+1

        s.sendto(b'RESET:20', (HOST, PORT))
        ready = select.select([s],[],[],1)

        if ready[0]:
            data,addr = s.recvfrom(1024)
            break
    print('received:',data)
    repeat = repeat+1
    total = total+i
    mean = total/repeat

print('Average:', mean)

s.close()