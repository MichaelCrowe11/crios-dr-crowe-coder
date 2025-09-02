"""
Authentication and Authorization System for CriOS Platform
JWT tokens, OAuth2, and role-based access control
"""

from datetime import datetime, timedelta
from typing import Optional, Dict, List
from fastapi import Depends, HTTPException, status, Security
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm, HTTPBearer, HTTPAuthorizationCredentials
from jose import JWTError, jwt
from passlib.context import CryptContext
from pydantic import BaseModel, EmailStr, Field
from sqlalchemy.orm import Session
import secrets
import stripe
from enum import Enum

# Configuration
SECRET_KEY = secrets.token_urlsafe(32)
ALGORITHM = "HS256"
ACCESS_TOKEN_EXPIRE_MINUTES = 30
REFRESH_TOKEN_EXPIRE_DAYS = 7

# Stripe configuration
stripe.api_key = "sk_live_your_stripe_key"  # Replace with actual key

# Password hashing
pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

# OAuth2 scheme
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/api/auth/token")
security = HTTPBearer()

# ============================================================================
# USER MODELS & TIERS
# ============================================================================

class UserTier(str, Enum):
    FREE = "free"
    RESEARCHER = "researcher"        # $99/month
    PROFESSIONAL = "professional"    # $499/month
    ENTERPRISE = "enterprise"        # $2,499/month
    ACADEMIC = "academic"           # $49/month (with verification)

class UserRole(str, Enum):
    USER = "user"
    RESEARCHER = "researcher"
    ADMIN = "admin"
    ENTERPRISE = "enterprise"

class TierLimits:
    """API rate limits and features per tier"""
    LIMITS = {
        UserTier.FREE: {
            "api_calls_per_day": 100,
            "agents_available": 10,
            "pipelines_per_month": 1,
            "compounds_per_run": 100,
            "storage_gb": 1,
            "support": "community",
            "features": ["basic_search", "molecule_viewer"]
        },
        UserTier.RESEARCHER: {
            "api_calls_per_day": 5000,
            "agents_available": 50,
            "pipelines_per_month": 10,
            "compounds_per_run": 10000,
            "storage_gb": 100,
            "support": "email",
            "features": ["all_free", "similarity_search", "clustering", "export"]
        },
        UserTier.PROFESSIONAL: {
            "api_calls_per_day": 50000,
            "agents_available": 194,
            "pipelines_per_month": 100,
            "compounds_per_run": 100000,
            "storage_gb": 1000,
            "support": "priority",
            "features": ["all_researcher", "custom_pipelines", "api_access", "collaboration"]
        },
        UserTier.ENTERPRISE: {
            "api_calls_per_day": -1,  # Unlimited
            "agents_available": 194,
            "pipelines_per_month": -1,  # Unlimited
            "compounds_per_run": -1,    # Unlimited
            "storage_gb": -1,           # Unlimited
            "support": "dedicated",
            "features": ["all_features", "white_label", "on_premise", "sla", "custom_agents"]
        },
        UserTier.ACADEMIC: {
            "api_calls_per_day": 10000,
            "agents_available": 194,
            "pipelines_per_month": 50,
            "compounds_per_run": 50000,
            "storage_gb": 500,
            "support": "email",
            "features": ["all_researcher", "academic_collaboration", "citation_exports"]
        }
    }

class User(BaseModel):
    id: str
    email: EmailStr
    username: str
    full_name: str
    tier: UserTier = UserTier.FREE
    role: UserRole = UserRole.USER
    is_active: bool = True
    created_at: datetime
    stripe_customer_id: Optional[str] = None
    stripe_subscription_id: Optional[str] = None
    api_key: Optional[str] = None
    organization: Optional[str] = None
    usage: Dict = Field(default_factory=dict)

class UserCreate(BaseModel):
    email: EmailStr
    username: str
    password: str
    full_name: str
    organization: Optional[str] = None

class UserLogin(BaseModel):
    username: str
    password: str

class Token(BaseModel):
    access_token: str
    refresh_token: str
    token_type: str = "bearer"
    tier: UserTier
    expires_in: int

# ============================================================================
# AUTHENTICATION FUNCTIONS
# ============================================================================

def verify_password(plain_password: str, hashed_password: str) -> bool:
    """Verify password against hash"""
    return pwd_context.verify(plain_password, hashed_password)

def get_password_hash(password: str) -> str:
    """Hash password"""
    return pwd_context.hash(password)

def create_access_token(data: dict, expires_delta: Optional[timedelta] = None) -> str:
    """Create JWT access token"""
    to_encode = data.copy()
    expire = datetime.utcnow() + (expires_delta or timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES))
    to_encode.update({"exp": expire, "type": "access"})
    return jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)

def create_refresh_token(data: dict) -> str:
    """Create JWT refresh token"""
    to_encode = data.copy()
    expire = datetime.utcnow() + timedelta(days=REFRESH_TOKEN_EXPIRE_DAYS)
    to_encode.update({"exp": expire, "type": "refresh"})
    return jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)

def create_api_key(user_id: str) -> str:
    """Create API key for programmatic access"""
    return f"crios_{'live' if True else 'test'}_{secrets.token_urlsafe(32)}"

async def get_current_user(credentials: HTTPAuthorizationCredentials = Security(security)) -> User:
    """Get current user from JWT token"""
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    
    try:
        token = credentials.credentials
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        user_id: str = payload.get("sub")
        tier: str = payload.get("tier", UserTier.FREE)
        
        if user_id is None:
            raise credentials_exception
            
        # TODO: Fetch user from database
        user = User(
            id=user_id,
            email=payload.get("email"),
            username=payload.get("username"),
            full_name=payload.get("full_name"),
            tier=tier,
            role=payload.get("role", UserRole.USER),
            created_at=datetime.now()
        )
        return user
        
    except JWTError:
        raise credentials_exception

def check_tier_limits(user: User, resource: str, amount: int = 1) -> bool:
    """Check if user has exceeded tier limits"""
    limits = TierLimits.LIMITS[user.tier]
    
    if resource == "api_calls":
        daily_calls = user.usage.get("api_calls_today", 0)
        limit = limits["api_calls_per_day"]
        if limit == -1:  # Unlimited
            return True
        return daily_calls + amount <= limit
    
    elif resource == "agents":
        return amount <= limits["agents_available"]
    
    elif resource == "pipelines":
        monthly_pipelines = user.usage.get("pipelines_this_month", 0)
        limit = limits["pipelines_per_month"]
        if limit == -1:  # Unlimited
            return True
        return monthly_pipelines + amount <= limit
    
    elif resource == "compounds":
        limit = limits["compounds_per_run"]
        if limit == -1:  # Unlimited
            return True
        return amount <= limit
    
    return False

def require_tier(minimum_tier: UserTier):
    """Decorator to require minimum tier for endpoint"""
    def decorator(user: User = Depends(get_current_user)):
        tier_order = {
            UserTier.FREE: 0,
            UserTier.ACADEMIC: 1,
            UserTier.RESEARCHER: 2,
            UserTier.PROFESSIONAL: 3,
            UserTier.ENTERPRISE: 4
        }
        
        if tier_order.get(user.tier, 0) < tier_order.get(minimum_tier, 0):
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail=f"This feature requires {minimum_tier} tier or higher. Current tier: {user.tier}"
            )
        return user
    return decorator

# ============================================================================
# STRIPE INTEGRATION
# ============================================================================

class StripeService:
    """Handle Stripe payments and subscriptions"""
    
    PRICE_IDS = {
        UserTier.RESEARCHER: "price_researcher_monthly",      # $99/month
        UserTier.PROFESSIONAL: "price_professional_monthly",  # $499/month
        UserTier.ENTERPRISE: "price_enterprise_monthly",      # $2,499/month
        UserTier.ACADEMIC: "price_academic_monthly",         # $49/month
    }
    
    @staticmethod
    async def create_customer(user: User) -> str:
        """Create Stripe customer"""
        customer = stripe.Customer.create(
            email=user.email,
            name=user.full_name,
            metadata={
                "user_id": user.id,
                "username": user.username,
                "organization": user.organization or ""
            }
        )
        return customer.id
    
    @staticmethod
    async def create_subscription(customer_id: str, tier: UserTier) -> Dict:
        """Create subscription for user"""
        price_id = StripeService.PRICE_IDS.get(tier)
        if not price_id:
            raise ValueError(f"No price configured for tier: {tier}")
        
        subscription = stripe.Subscription.create(
            customer=customer_id,
            items=[{"price": price_id}],
            trial_period_days=14,  # 14-day free trial
            metadata={"tier": tier}
        )
        return {
            "subscription_id": subscription.id,
            "status": subscription.status,
            "current_period_end": subscription.current_period_end
        }
    
    @staticmethod
    async def create_checkout_session(user: User, tier: UserTier, success_url: str, cancel_url: str) -> str:
        """Create Stripe checkout session"""
        price_id = StripeService.PRICE_IDS.get(tier)
        
        session = stripe.checkout.Session.create(
            payment_method_types=["card"],
            line_items=[{
                "price": price_id,
                "quantity": 1,
            }],
            mode="subscription",
            success_url=success_url,
            cancel_url=cancel_url,
            customer_email=user.email,
            metadata={
                "user_id": user.id,
                "tier": tier
            },
            subscription_data={
                "trial_period_days": 14,
                "metadata": {"tier": tier}
            }
        )
        return session.url
    
    @staticmethod
    async def cancel_subscription(subscription_id: str) -> bool:
        """Cancel subscription"""
        try:
            stripe.Subscription.delete(subscription_id)
            return True
        except Exception as e:
            print(f"Error canceling subscription: {e}")
            return False
    
    @staticmethod
    async def create_usage_record(subscription_item_id: str, quantity: int) -> Dict:
        """Record usage for metered billing"""
        usage_record = stripe.SubscriptionItem.create_usage_record(
            subscription_item_id,
            quantity=quantity,
            timestamp=int(datetime.now().timestamp())
        )
        return usage_record

# ============================================================================
# OAUTH2 PROVIDERS
# ============================================================================

class OAuth2Service:
    """Handle OAuth2 authentication with external providers"""
    
    @staticmethod
    async def google_auth(token: str) -> Dict:
        """Authenticate with Google OAuth2"""
        # Verify Google token and get user info
        # Implementation would use google-auth library
        pass
    
    @staticmethod
    async def github_auth(token: str) -> Dict:
        """Authenticate with GitHub OAuth2"""
        # Verify GitHub token and get user info
        # Implementation would use GitHub API
        pass
    
    @staticmethod
    async def microsoft_auth(token: str) -> Dict:
        """Authenticate with Microsoft OAuth2"""
        # Verify Microsoft token and get user info
        # Implementation would use MSAL library
        pass

# ============================================================================
# API KEY MANAGEMENT
# ============================================================================

class APIKeyService:
    """Manage API keys for programmatic access"""
    
    @staticmethod
    async def create_api_key(user: User, name: str = "Default") -> Dict:
        """Create new API key for user"""
        api_key = create_api_key(user.id)
        
        # Store in database with user association
        # TODO: Implement database storage
        
        return {
            "api_key": api_key,
            "name": name,
            "created_at": datetime.now().isoformat(),
            "tier": user.tier,
            "rate_limit": TierLimits.LIMITS[user.tier]["api_calls_per_day"]
        }
    
    @staticmethod
    async def validate_api_key(api_key: str) -> Optional[User]:
        """Validate API key and return associated user"""
        # TODO: Implement database lookup
        # For now, return mock user
        if api_key.startswith("crios_"):
            return User(
                id="api_user_123",
                email="api@example.com",
                username="api_user",
                full_name="API User",
                tier=UserTier.PROFESSIONAL,
                role=UserRole.USER,
                created_at=datetime.now(),
                api_key=api_key
            )
        return None
    
    @staticmethod
    async def revoke_api_key(api_key: str) -> bool:
        """Revoke API key"""
        # TODO: Implement database update
        return True

# ============================================================================
# USAGE TRACKING
# ============================================================================

class UsageTracker:
    """Track and limit API usage per tier"""
    
    @staticmethod
    async def track_api_call(user: User, endpoint: str, tokens_used: int = 0):
        """Track API call for billing and rate limiting"""
        # TODO: Implement Redis-based tracking
        pass
    
    @staticmethod
    async def track_pipeline_run(user: User, pipeline: str, compounds: int):
        """Track pipeline execution for billing"""
        # TODO: Implement database storage
        pass
    
    @staticmethod
    async def get_usage_report(user: User, period: str = "month") -> Dict:
        """Get usage report for user"""
        return {
            "period": period,
            "api_calls": 1234,
            "pipelines_run": 5,
            "compounds_generated": 50000,
            "storage_used_gb": 12.5,
            "cost_estimate": 99.00
        }

# Example: Protected endpoint with tier checking
async def protected_endpoint(
    user: User = Depends(require_tier(UserTier.RESEARCHER))
):
    """Example endpoint requiring Researcher tier or higher"""
    return {"message": f"Hello {user.full_name}, you have {user.tier} access"}