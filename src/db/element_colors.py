import sqlite3
from logging import getLogger

from constants import NAME


class ElementColors:
    def __init__(self, connection):
        self.logger = getLogger(NAME)
        self.connection = connection
        self.cursor = connection.cursor()
        self.init()

        self.logger.info("ElementColors initialized.")

    def init(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS element_colors (
                server_id TEXT,
                element TEXT,
                color TEXT,
                PRIMARY KEY (server_id, element)
            )
        ''')
        self.connection.commit()

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
        for row in results:
            atom_palette[int(row["element"])] = tuple(int(rgb) for rgb in row["color"].split(','))

        return atom_palette
    
    def clear(self, server_id):
        self.cursor.execute('''
        DELETE FROM element_colors WHERE server_id = ?;
        ''', (server_id,))
        self.connection.commit()