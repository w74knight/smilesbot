import sqlite3

class DatabaseHandler:
    def __init__(self, db_name='database.db'):
        self.connection = sqlite3.connect(db_name)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()
        self.create_table()

    def create_table(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS server_settings (
                server_id TEXT PRIMARY KEY,
                prefix TEXT,
                auto_smile BOOLEAN
            )
        ''')
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS element_colors (
                server_id TEXT,
                element TEXT,
                color TEXT,
                PRIMARY KEY (server_id, element)
            )
        ''')
        self.connection.commit()

    def set_server_setting(self, server_id, prefix):
        self.cursor.execute('''
            INSERT INTO server_settings (server_id, prefix)
            VALUES (?, ?)
            ON CONFLICT(server_id) DO UPDATE SET
            prefix=excluded.prefix
        ''', (server_id, prefix))
        self.connection.commit()

    def update_server_setting(self, server_id, key, value):
        self.cursor.execute(f'''
            UPDATE server_settings
            SET {key}=?
            WHERE server_id=?
        ''', (value, server_id))
        self.connection.commit()

    def get_server_setting(self, server_id):
        self.cursor.execute('''
            SELECT * FROM server_settings WHERE server_id=?
        ''', (server_id,))
        result = self.cursor.fetchone()
        return dict(result) if result else {}
    
    def set_element_color(self, server_id, element, color):
        self.cursor.execute('''
            INSERT INTO element_colors (server_id, element, color)
            VALUES (?, ?, ?)
            ON CONFLICT(server_id, element) DO UPDATE SET
            color=excluded.color
        ''', (server_id, element, color))
        self.connection.commit()

    def get_element_colors(self, server_id):
        self.cursor.execute('''
            SELECT element, color FROM element_colors WHERE server_id=?
        ''', (server_id,))
        results = self.cursor.fetchall()
        atom_palette = {}
        if results:
            for row in results:
                atom_palette[int(row["element"])] = tuple(int(rgb) / 255 for rgb in row["color"].split(', '))
        return atom_palette
    
    def create_server_table(self, server_id):
        self.cursor.execute('''
            INSERT OR IGNORE INTO server_settings (server_id, prefix, auto_smile)
            VALUES (?, ?, ?)
        ''', (server_id, '/', False))
        self.connection.commit()

    def clear_color(self, server_id):
        self.cursor.execute('''
        DELETE FROM element_colors WHERE server_id = ?;
        ''', (server_id,))
        self.connection.commit()

    def clear_settings(self, server_id):
        self.cursor.execute('''
        DELETE FROM server_settings WHERE server_id = ?;
        ''', (server_id,))
        self.connection.commit()

    def close(self):
        self.connection.close()