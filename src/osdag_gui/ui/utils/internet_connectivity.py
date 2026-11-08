import socket


class InternetConnectivity:
    def is_online(self, timeout=3):
        try:
            socket.create_connection(("1.1.1.1", 53), timeout=timeout).close()
            return True
        except OSError:
            return False
