import sqlite3
from logging import getLogger
import json
from constants import NAME
from constants import DISCORD_DARK_JSON as discord_dark

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

    def set_element_defaults(self, server_id, DISCORD_DARK_JSON):
        for element, color in DISCORD_DARK_JSON.items():
            color_json = json.dumps(color)
            self.cursor.execute('''
                INSERT OR REPLACE INTO element_colors (server_id, element, color)
                VALUES (? , ?, ?)
            ''', (server_id, str(element), color_json))
        self.connection.commit()

    def get_element_colors(self, server_id):
        self.cursor.execute('''
            SELECT element, color FROM element_colors WHERE server_id=?
        ''', (server_id,))
        results = self.cursor.fetchall()
        atom_palette = {}
        for row in results:
            element = int(row["element"])
            color_list = json.loads(row["color"])  # parse JSON back to list
            atom_palette[element] = tuple(color_list)  # convert list to tuple

        return atom_palette
    
    def clear(self, server_id):
        self.cursor.execute('''
        DELETE FROM element_colors WHERE server_id = ?;
        ''', (server_id,))
        self.connection.commit()