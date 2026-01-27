
/*

## USER

Is there a simple, way how to implement simple fultext search ? 
(resp similarity search, if I do not know exactly how the substance is labeled, and I want to lost all with similar names? 
What are the possibilites? I guess there are libraries and web services for that? mention them. 
But if we can implement something simple oursefl it would be good. 
Notice that it should be only on very short strings. (similarty seach, you see I do a lot of typos, it should be typo tolerant)

---

## Gemini

For short strings like chemical names and formulas, the most robust and "classic" way to handle typos and similarity is the Levenshtein Distance algorithm. It calculates the minimum number of single-character edits (insertions, deletions, substitutions) required to change one word into another.
Since you want to keep it simple and local, I will show you:
A pure JavaScript implementation of a similarity search.
How to integrate it into your Compound search logic.
A mention of libraries if you want to go "pro" later.
1. The Similarity Implementation (Searcher.js)
Here is a simple class you can use. It normalizes strings (lowercase) and calculates a similarity score between 0 and 1.


### 3. What are the possibilities?

#### A. DIY Approaches (What we just did)
*   **Levenshtein (Edit Distance):** Best for typos (e.g., "Sulfur" vs "Sulfer").
*   **N-Grams:** Break words into chunks (e.g., "Acid" -> `ac`, `ci`, `id`). If two words share many chunks, they are similar. Better for longer strings.
*   **Phonetic Algorithms (Soundex/Metaphone):** Matches words that *sound* the same. Great for names (e.g., "Smith" vs "Smyth"), but less useful for chemical formulas like "H2SO4".

#### B. JavaScript Libraries (Recommended for Browser integration later)
*   **[Fuse.js](https://fusejs.io/):** The gold standard for lightweight fuzzy search in JS. It requires zero dependencies and works great with arrays of objects. It uses a more sophisticated version of the logic I wrote above.
*   **[FlexSearch](https://github.com/nextapps-de/flexsearch):** If your database grows to thousands of compounds, this is the fastest library available.

#### C. Web Services (For "Big Data")
*   **PubChem Autocomplete API:** You can send a partial string to PubChem, and they return the "official" chemical name and CID.
    *   *URL Example:* `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/glucose/description/JSON`
*   **Algolia:** A professional search-as-a-service. Expensive, but offers "google-like" search speed and typo tolerance.

### My Recommendation for you:
Since you are working with a limited list of compounds (hundreds or a few thousand), **Fuse.js** is your best bet for the web version. 

For the **Node.js terminal version**, the `ChemicalSearcher` class I provided above is perfectly sufficient, adds zero dependencies, and is very "typo-tolerant" for exactly the kind of mistakes you mentioned. 

**Try running this:**
Search for `"Yperite"`—it should find `"Sulfur mustard"` because it's in the `notes` column, or search for `"H2SO3"`—it should suggest `"H2SO4"` with a high similarity score.

*/

export class StringSearcher {
    /**
     * Calculates the Levenshtein distance between two strings.
     */
    static levenshtein(a, b) {
        const matrix = [];
        for (let i = 0; i <= b.length; i++) matrix[i] = [i];
        for (let j = 0; j <= a.length; j++) matrix[0][j] = j;

        for (let i = 1; i <= b.length; i++) {
            for (let j = 1; j <= a.length; j++) {
                if (b.charAt(i - 1) === a.charAt(j - 1)) {
                    matrix[i][j] = matrix[i - 1][j - 1];
                } else {
                    matrix[i][j] = Math.min(
                        matrix[i - 1][j - 1] + 1, // substitution
                        matrix[i][j - 1] + 1,     // insertion
                        matrix[i - 1][j] + 1      // deletion
                    );
                }
            }
        }
        return matrix[b.length][a.length];
    }

    /**
     * Returns a score from 0 to 1 (1 being exact match)
     */
    static compare(str1, str2) {
        const s1 = str1.toLowerCase();
        const s2 = str2.toLowerCase();
        if (s1 === s2) return 1.0;
        
        const distance = this.levenshtein(s1, s2);
        const maxLength = Math.max(s1.length, s2.length);
        return 1.0 - distance / maxLength;
    }

    /**
     * Searches a collection for a query string.
     * @param {Array} collection - Array of objects
     * @param {string} query - The user input
     * @param {Array} keys - Which object keys to search (e.g. ['name', 'formula'])
     * @param {number} threshold - Minimum similarity (0.4 - 0.6 is usually good)
     */
    static search(collection, query, keys, threshold = 0.5) {
        return collection
            .map(item => {
                // Find the best match among all specified keys
                const bestScore = Math.max(...keys.map(key => {
                    const val = String(item[key] || "");
                    // Check if it contains the substring exactly (bonus)
                    if (val.toLowerCase().includes(query.toLowerCase())) return 0.9;
                    return this.compare(query, val);
                }));
                return { item, score: bestScore };
            })
            .filter(result => result.score >= threshold)
            .sort((a, b) => b.score - a.score); // Highest score first
    }
}

export default StringSearcher;