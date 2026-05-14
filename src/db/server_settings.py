from logging import getLogger
from sqlite3 import Connection

from constants import NAME

DEFAULT_SERVER_SETTINGS = {
    "prefix": "/",
    "auto_smile": False
}

class ServerSettings:
    def __init__(self, connection):
        self.logger = getLogger(NAME)

        self.connection = connection
        self.cursor = connection.cursor()
        self.init()

        self.logger.info("ServerSettings initialized.")

    def init(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS server_settings (
                server_id TEXT PRIMARY KEY,
                prefix TEXT,
                auto_smile BOOLEAN
            )
        ''')
        self.connection.commit()

    def set_server_setting(self, server_id, key, value):
        query = f'''
            INSERT INTO server_settings (server_id, {key})
            VALUES (?, ?)
            ON CONFLICT(server_id) DO UPDATE SET
            {key}=excluded.{key}
        '''

        # Execute the query with the server_id and value
        self.cursor.execute(query, (server_id, value))
        self.connection.commit()

    def get_server_setting(self, server_id):
        self.cursor.execute('''
            SELECT * FROM server_settings WHERE server_id=?
        ''', (server_id,))
        result = self.cursor.fetchone()
        return dict(result) if result else {}
    
    def clear(self, server_id):
        self.cursor.execute('''
        DELETE FROM server_settings WHERE server_id = ?;
        ''', (server_id,))
        self.connection.commit()

    def set_server_defaults(self, server_id):
        self.cursor.execute('''
            INSERT INTO server_settings (server_id, prefix, auto_smile)
            VALUES (?, ?, ?)
            ON CONFLICT(server_id) DO UPDATE SET
                prefix=excluded.prefix,
                auto_smile=excluded.auto_smile
        ''', (
            server_id,
            DEFAULT_SERVER_SETTINGS["prefix"],
            DEFAULT_SERVER_SETTINGS["auto_smile"]
        ))
        self.connection.commit()