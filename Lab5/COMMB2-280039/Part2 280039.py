import socket

HOST = "128.178.151.64"  # "tcpip.epfl.ch"
PORT = 5003

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.connect((HOST, PORT))

# try:
# finally:
# sock.close()

i=0
message = "CMD_floodme"
# message = "CMD_short:0"
print('sending:', message)
sock.sendall(message.encode())

while True:
    data = sock.recv(10000).decode()
    # data = sock.recv(18).decode()
    i = i+1
    if data:
        print('received:', data)
    else:
        print('No more data from the server')
        break

print('recv invoked times:',i)

sock.close()
