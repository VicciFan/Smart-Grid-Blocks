import socket

#HOST = "localhost"
HOST = "128.178.151.64" #"tcpip.epfl.ch"
PORT = 5002   

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.connect((HOST,PORT))

#try:
#finally:
#	sock.close()

message = "O Romeo, Romeo! wherefore art thou Romeo?"
print("sending:", message)
sock.sendall(message.encode())
	
received = 0
expected = len(message)
	
while received < expected:
	data = sock.recv(16).decode()
	received += len(data)
	print ('received:', data)

while True:
	pass
