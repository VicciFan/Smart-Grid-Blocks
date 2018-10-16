import socket

HOST = 'localhost'           
PORT = 5002  

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
sock.bind((HOST, PORT))
sock.listen(1)
while True:
	connection, addr = sock.accept()
	while True:
		data = connection.recv(16).decode()
		print("received:", data)
		if data:
			connection.sendall(data.encode())
		else:
			print("No more data from", addr)
			break
	
	connection.close()
