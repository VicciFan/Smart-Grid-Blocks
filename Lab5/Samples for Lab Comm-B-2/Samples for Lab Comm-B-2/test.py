import socket

import time

HOST = '2002:80b2:9c51::ffff:901'  # The remote host

PORT = 5004  # The same port as used by the server

s = socket.socket(socket.AF_INET6, socket.SOCK_DGRAM)

iter = 0

sum = 0

while iter <= 60:

iter = iter + 1

i = 0

while True:

i = i + 1

s.sendto(b'RESET:20', (HOST, PORT))

data, addr = s.recvfrom(1024)

if data:

break

else:

time.sleep(1)

print('data received at', i)

print(data)

sum = sum + i
avg = sum / 60

print('average received time:', avg)

s.close()