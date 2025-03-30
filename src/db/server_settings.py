class ServerSettings:
    def __init__(self, connection):
        self.connection = connection
        self.cursor = connection.cursor()
        self.init()

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