import sqlite3
from logging import getLogger

from constants import NAME, SMILE_BG

DEFAULT_RENDER_OPTIONS = {
    "background_color": ",".join(str(c) for c in SMILE_BG),
    "includeAtomNumbers": False,
    "addStereoAnnotations": False,
    "explicitMethyl": False,
    "atomLabelDeuteriumTritium": False,
    "dummiesAreAttachments": False
}

class RenderOptions:
    def __init__(self, connection):
        self.logger = getLogger(NAME)
        self.connection = connection
        self.cursor = connection.cursor()
        self.init()

    def init(self):
        self.logger = getLogger(NAME)
        self.init()

        self.logger.info("RenderOptions initialized.")

    def init(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS render_options (
                server_id TEXT PRIMARY KEY,
                background_color TEXT,
                includeAtomNumbers BOOLEAN,
                addStereoAnnotations BOOLEAN,
                explicitMethyl BOOLEAN,
                atomLabelDeuteriumTritium BOOLEAN,
                dummiesAreAttachments BOOLEAN
            )
        ''')
        self.connection.commit()

    def get_render_option_keys(self) -> set[str]:
        self.cursor.execute("PRAGMA table_info(render_options)")
        rows = self.cursor.fetchall()
        # Exclude non-rendering columns
        skip = {"server_id", "background_color"}
        return {row["name"] for row in rows if row["name"] not in skip}

    # Get render options
    def get_render_option(self, server_id):
        self.cursor.execute('''
            SELECT includeAtomNumbers, addStereoAnnotations, explicitMethyl, 
                   atomLabelDeuteriumTritium, dummiesAreAttachments
            FROM render_options WHERE server_id=?
        ''', (server_id,))
        result = self.cursor.fetchone()
        return dict(result) if result else DEFAULT_RENDER_OPTIONS

    # Setter/Getter for background color
    def set_bgcolor(self, server_id, color):
        self.cursor.execute('''
            INSERT INTO render_options (server_id, background_color)
            VALUES (?, ?)
            ON CONFLICT(server_id) DO UPDATE SET
            background_color=excluded.background_color
        ''', (server_id, color))
        self.connection.commit()

    def get_bgcolor(self, server_id):
        self.cursor.execute('''
            SELECT background_color FROM render_options WHERE server_id=?
        ''', (server_id,))
        result = self.cursor.fetchone()
        
        if result and result["background_color"]:
            bg_color = result["background_color"].split(",")
            return tuple(int(c) for c in bg_color)

        return SMILE_BG

    def set(self, server_id, key, value):
        query = f'''
            INSERT INTO render_options (server_id, {key})
            VALUES (?, ?)
            ON CONFLICT(server_id) DO UPDATE SET
            {key}=excluded.{key}
        '''
        self.cursor.execute(query, (server_id, value))
        self.connection.commit()

    def get(self, server_id, key):
        self.cursor.execute(f'''SELECT {key} FROM render_options WHERE server_id=?''', (server_id,))
        result = self.cursor.fetchone()
        return result[key] if result else False
    
    def clear(self, server_id):
        self.cursor.execute('''
        DELETE FROM render_options WHERE server_id = ?;
        ''', (server_id,))
        self.connection.commit()

    def set_render_defaults(self, server_id):
        self.cursor.execute('''
            INSERT INTO render_options (server_id, background_color, includeAtomNumbers, addStereoAnnotations, 
                                       explicitMethyl, atomLabelDeuteriumTritium, dummiesAreAttachments)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(server_id) DO UPDATE SET
                background_color=excluded.background_color,
                includeAtomNumbers=excluded.includeAtomNumbers,
                addStereoAnnotations=excluded.addStereoAnnotations,
                explicitMethyl=excluded.explicitMethyl,
                atomLabelDeuteriumTritium=excluded.atomLabelDeuteriumTritium,
                dummiesAreAttachments=excluded.dummiesAreAttachments
        ''', (
            server_id,
            DEFAULT_RENDER_OPTIONS["background_color"],
            DEFAULT_RENDER_OPTIONS["includeAtomNumbers"],
            DEFAULT_RENDER_OPTIONS["addStereoAnnotations"],
            DEFAULT_RENDER_OPTIONS["explicitMethyl"],
            DEFAULT_RENDER_OPTIONS["atomLabelDeuteriumTritium"],
            DEFAULT_RENDER_OPTIONS["dummiesAreAttachments"]
        ))
        self.connection.commit()