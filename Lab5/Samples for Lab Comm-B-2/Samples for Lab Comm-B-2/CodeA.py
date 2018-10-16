import socket
HOST = 'localhost' # The remote host
PORT = 50007 # The same port as used by the server
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.sendto(b'Hello, Romeo ', (HOST,PORT) )
s.close()
print('Message sent')
