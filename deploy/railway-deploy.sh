#!/bin/bash

# CriOS Discovery Engine - Railway Deployment Script
# Quick deployment helper for Railway.app

set -e

echo "🚂 CriOS Discovery Engine - Railway Deployment Helper"
echo "=================================================="
echo ""

# Check if Railway CLI is installed
if ! command -v railway &> /dev/null; then
    echo "❌ Railway CLI not found. Installing..."
    echo ""
    echo "Please run one of these commands first:"
    echo "  npm i -g @railway/cli"
    echo "  or"
    echo "  bash <(curl -fsSL cli.new)"
    echo ""
    exit 1
fi

echo "✅ Railway CLI found"
echo ""

# Check if logged in
if ! railway whoami &> /dev/null; then
    echo "🔐 Please login to Railway..."
    railway login
fi

echo "✅ Logged in to Railway"
echo ""

# Create new project or link existing
read -p "Create new project or link existing? (new/link): " choice

if [ "$choice" = "new" ]; then
    echo "📦 Creating new Railway project..."
    railway init
elif [ "$choice" = "link" ]; then
    echo "🔗 Linking to existing project..."
    railway link
else
    echo "❌ Invalid choice"
    exit 1
fi

echo ""
echo "🔧 Setting up services..."
echo ""

# Function to add a service
add_service() {
    local service_name=$1
    echo "Adding $service_name..."
    railway service create $service_name || echo "Service already exists"
}

# Add PostgreSQL
echo "📊 Adding PostgreSQL database..."
railway add --database postgresql || echo "PostgreSQL already added"
echo ""

# Add Redis
echo "💾 Adding Redis cache..."
railway add --database redis || echo "Redis already added"
echo ""

# Environment variables setup
echo "🔐 Setting up environment variables..."
echo ""

read -p "Enter your Anthropic API key: " anthropic_key

if [ -z "$anthropic_key" ]; then
    echo "⚠️  No API key provided. You'll need to set ANTHROPIC_API_KEY later."
else
    railway variables set ANTHROPIC_API_KEY="$anthropic_key"
    echo "✅ Anthropic API key set"
fi

echo ""
echo "🚀 Deploying services..."
echo ""

# Deploy backend
echo "📦 Deploying Backend API..."
railway up -d --service backend || railway up -d

echo ""
echo "📦 Deploying Celery Worker..."
railway up -d --service celery-worker || echo "Celery service deployment skipped"

echo ""
echo "📦 Deploying Frontend..."
railway up -d --service frontend || echo "Frontend deployment skipped"

echo ""
echo "✅ Deployment initiated!"
echo ""
echo "📊 View deployment status:"
echo "   railway status"
echo ""
echo "📝 View logs:"
echo "   railway logs"
echo ""
echo "🌐 Open project dashboard:"
echo "   railway open"
echo ""
echo "📚 For detailed setup, see RAILWAY_DEPLOYMENT.md"
echo ""
echo "🎉 Done! Your CriOS Discovery Engine is deploying to Railway!"
