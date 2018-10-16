import websocket
import socket

from websocket import create_connection

ws = create_connection("ws://tcpip.epfl.ch:5006")
ws.send("CMD_short:0")
# ws.send("CMD floodme")

while True:
    data = ws.recv_data(10000)
    if data:
        print('received:', data)
    else:
        break

ws.close()