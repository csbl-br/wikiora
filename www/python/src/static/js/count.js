function fetchListsEnriched() {
  return fetch('/api/lists_enriched')
      .then(response => response.json())
      .catch(error => {
          console.error('Error fetching lists enriched:', error);
          return { lists_enriched: 0 };
      });
}

async function updateCounter() {
  const data = await fetchListsEnriched();
  const newCount = data.lists_enriched;
  smoothIncrement(0, newCount);
}

function smoothIncrement(start, end) {
  const duration = 1000; // Duration of the animation in milliseconds
  const increment = (end - start) / (duration / 10);
  let current = start;
  const interval = setInterval(() => {
      current += increment;
      document.getElementById('lists-enriched-counter').textContent = Math.floor(current);
      if (current >= end) {
          clearInterval(interval);
          document.getElementById('lists-enriched-counter').textContent = end; // Ensure final value is set accurately
      }
  }, 10);
}


// Call updateCounter initially to set the counter value
updateCounter();

// Optionally, you can set an interval to update the counter periodically
// setInterval(updateCounter, 60000); // Update every minute
